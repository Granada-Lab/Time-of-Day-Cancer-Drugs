function ToD_growth_plots_parameter_extraction(date,cellline,drug,timepoints,treat1,treat2,finaltimepoint,intervals,v,color)

%Carolin Ector, 23.08.2023
%Time-of-Day-Cancer-Drugs Manuscript Fig.4b,c

%%Function reads and normalizes time-series growth data from Time-of-Day treatment experiments stored in excel sheets
%%Generation of growth plots and Time-of-Day (ToD) response curves

%input: stored in '[date]_ToD_workspace.mat'
% date: date of the experiment being analyzed
% cellline: names of the cell lines being analysed
% drug: names of drugs used for time-of-day treatments
% treat1: time of treatment 1 (h post start of recording)
% treat2: time of treatment 2 (h post start of recording)
% finaltimepoint: timepoint used for evaluation of final responses (h post treatment)
% intervals: imaging time intervals (h)
% v: initial column for first treatment time point for each drug. example: v(:,1) = column of initial treatment timepont for drug{1}
% color: color for different time points

%Define variables
channel = {'Conf';'CellNr'}; %confluency = brightfield channel, cell number = red fluorescent channel from live-imaging
yaxisnames = {'Confluency';'Cell Number'};
cc = numel(cellline);
dd = numel(drug);
ii = numel(timepoints);

%% Load and normalize data
experiment = str2double(date);

for a = 1:2 %loop a channel

    for c = 1:cc %loop c cell lines

        % Load excel file where data is stored
        pathtoexcel = append('ToD_experiments_raw_data/',date,'_ToD_',cellline{c},'.xlsx');
        sheet1 = append(channel{a},'_1');
        sheet2 = append(channel{a},'_2');
        [num1] = xlsread(pathtoexcel,sheet1); %data from plate 1
        [num2] = xlsread(pathtoexcel,sheet2); %data from plate 2

        %account for potentially different recording lengths across the two assayed plates
        length_array1 = size(num1,1); length_array2 = size(num2,1);
        if length_array1 > length_array2
            num1 = num1(1:length_array2,:);
        elseif length_array1 < length_array2
            num2 = num2(1:length_array1,:);
        end

        for d = 1:dd %loop drug

            %define order how to load columns (data organization is partially different for the different experiments
            %goal: load columns serially by timepoint (24 to 48h)
            if experiment == 20221107
                if d == 8 || d == 9
                    control_col = (51:3:66); %initial column for control
                    control_colorder = [control_col(3),control_col(2),control_col(1),control_col(6),control_col(5),control_col(4)];
                    if d == 8 %adavosertib
                        treat_colorder = control_colorder+1;
                    elseif d == 9 %cisplatin (separate NaCl control)
                        treat_colorder = control_colorder+2;
                    end
                else
                    v = (4:10); %initial column respective drug
                    treat_col = (v(d):8:50);  %initial column for drug-treated conditions
                    treat_colorder =  [treat_col(3),treat_col(2),treat_col(1),treat_col(6),treat_col(5),treat_col(4)];
                    control_col = (3:8:43);
                    control_colorder = [control_col(3),control_col(2),control_col(1),control_col(6),control_col(5),control_col(4)];
                end
            elseif experiment == 20220209
                treat_colorder = [v(d),v(d)-8,v(d)-16,v(d)+8,v(d)+16,v(d)+24];
                control_colorder = [19,11,3,27,35,43];
            elseif experiment == 20220215 || experiment == 20220228
                treat_colorder = (v(d):8:50);
                control_colorder=(3:8:43);
            elseif experiment == 20211001 || experiment == 20211027
                treat_colorder = (v:2:18);
                control_colorder = (v-1:2:17);
                if experiment == 20211001
                    num1 = num1(1:134,:);
                    num2 = num2(1:134,:);
                elseif experiment == 20211027
                    num1 = num1(1:104,:);
                    num2 = num2(1:104,:);
                end
            end

            for i = 1:ii %loop i timepoints

                for r = 1:2 %loop r plate replicates

                    if r == 1
                        num = num1;
                    elseif r == 2
                        num = num2;
                    end

                    %load x-values (=time)
                    xdata(:,r)=rmmissing(num(:,2));

                    %load y-data timepoints
                    yval_raw(:,r)=rmmissing(num(:,treat_colorder(i)));
                    ctrl_raw(:,r)=rmmissing(num(:,control_colorder(i)));

                    %smooth y-data
                    yval_smooth(:,r) = smoothdata(yval_raw(:,r));
                    ctrl_smooth(:,r) = smoothdata(ctrl_raw(:,r));

                    %other smoothing option
                    %                     yval_smooth(:,r) = smooth(yval_raw(:,r),0.35,'rloess');
                    %                     ctrl_smooth(:,r) = smooth(ctrl_raw(:,r),0.35,'rloess');

                    % 1. normalize ydata to timepoint 0
                    yval_normx0(:,r) = yval_smooth(:,r)/yval_smooth(1,r);
                    ctrlval_normx0(:,r) = ctrl_smooth(:,r)/ctrl_smooth(1,r);

                    % 2. normalize ydata to time of treatment (=ToT)
                    %find Time of Treatment (ToT) in x-data
                    time = xdata(:,r);
                    data = yval_smooth(:,r);
                    control = ctrl_smooth(:,r);

                    if experiment == 20220228 && c ~= 1
                        if treat1 == 51
                            treat1 = treat1 - 1.5;
                        end
                    end

                    if(i>=1 && i<=3)
                        ToT = time((time >= (treat1-1)) & (time <= (treat1+2)));
                    else
                        ToT = time((time >= (treat2-1)) & (time <= (treat2+2)));
                    end

                    % find y-value at x=ToT to normalize data to that point
                    ToT = min(ToT);
                    yval_ToT = data(time==ToT);

                    % find corresponding row numbers to shorten the dataset (exclude timepoints before ToT)
                    row_start = find(data == yval_ToT);
                    row_end = row_start + finaltimepoint/intervals;

                    %shorten datasets
                    ydata_ToT = data(row_start:row_end);
                    ctrl_ToT = control(row_start:row_end);
                    time_ToT = time(row_start:row_end);

                    %normalize shorter datasets to y-values at ToT
                    yval_norm_ToT(:,r) = ydata_ToT/ydata_ToT(time_ToT == time_ToT(1,:));
                    ctrl_norm_ToT(:,r) = ctrl_ToT/ctrl_ToT(time_ToT == time_ToT(1,:));

                    % 3. normalize ydata_normToT to control over time
                    yval_normctrl(:,r)=yval_norm_ToT(:,r)./ctrl_norm_ToT(:,r);

                    % 4. extract final responses for response curve
                    yval_finalTP(i,r) = median(yval_norm_ToT(finaltimepoint/intervals-3:end,r));
                    ctrlval_finalTP(i,r) = median(ctrl_norm_ToT(finaltimepoint/intervals-3:end,r));

                    %normalize all timepoints to TP0
                    yvaluesnormtoTP0(i,r) = yval_finalTP(i,r)./yval_finalTP(1,r);
                    ctrlvaluesnormtoTP0(i,r) = ctrlval_finalTP(i,r)./ctrlval_finalTP(1,r);

                    % 1.2 calculate maximal growth inhibition
                    inhibition(i,r) = mean(yval_normctrl((end-3:end),r));

                    valuestoclear1 = {'time';'data';'control';'ToT';'yval_ToT';
                        'row_start';'row_end';'ydata_ToT';'ctrl_ToT';'time_ToT'};
                    clear(valuestoclear1{:});

                end

                %take mean values of duplicates and save values for figures
                ydata_raw(:,i) = mean(yval_smooth, 2); %figure1
                ydata_raw_std(:,i) = std(yval_smooth,[],2);

                ydata_normx0(:,i) = mean(yval_normx0, 2); %figure2
                ydata_normx0_std(:,i) = std(yval_normx0,[],2);

                response_finalTP(i,d) = mean(yval_finalTP(i,:), 2); %figure 8
                response_finalTP_std(i,d) = std(yval_finalTP(i,:),[],2);

                if experiment == 20221107 && d == 8
                    ctrl_raw(:,i) = ctrl_smooth(:,1);
                    ctrl_raw_std(:,i) = (0);
                    ctrl_normx0(:,i) = ctrlval_normx0(:,1);
                    ctrl_normx0_std(:,i) = (0);
                    ctrl_finalTP(i,d) = ctrlval_finalTP(i,1);
                    ctrl_finalTP_std(i,d) = (0);
                    ctrl_normtoTP0(i,d) = ctrlvaluesnormtoTP0(i,1);
                    ctrl_normtoTP0_std(i,d) = (0);
                elseif experiment == 20221107 && d == 9
                    ctrl_raw(:,i) = ctrl_smooth(:,2);
                    ctrl_raw_std(:,i) = (0);
                    ctrl_normx0(:,i) = ctrlval_normx0(:,2);
                    ctrl_normx0_std(:,i) = (0);
                    ctrl_finalTP(i,d) = ctrlval_finalTP(i,2);
                    ctrl_finalTP_std(i,d) = (0);
                    ctrl_normtoTP0(i,d) = ctrlvaluesnormtoTP0(i,2);
                    ctrl_normtoTP0_std(i,d) = (0);
                else
                    ctrl_raw(:,i) = mean(ctrl_smooth, 2);
                    ctrl_raw_std(:,i) = std(ctrl_smooth,[],2);
                    ctrl_normx0(:,i) = mean(ctrlval_normx0, 2);
                    ctrl_normx0_std(:,i) = std(ctrlval_normx0,[],2);
                    ctrl_finalTP(i,d) = mean(ctrlval_finalTP(i,:),2); %figure 8
                    ctrl_finalTP_std(i,d) = std(ctrlval_finalTP(i,:),[],2);
                    ctrl_normtoTP0(i,d) = mean(ctrlvaluesnormtoTP0(i,:),2); %figure0.1
                    ctrl_normtoTP0_std(i,d) = std(ctrlvaluesnormtoTP0(i,:),[],2);
                end

                ydata_normToT(:,i) = mean(yval_norm_ToT,2); %figure3
                ydata_normToT_std(:,i) = std(yval_norm_ToT,[],2);
                ctrl_normToT(:,i) = mean(ctrl_norm_ToT,2);
                ctrl_normToT_std(:,i) = std(ctrl_norm_ToT,[],2);

                y_normctrl(:,i) = mean(yval_normctrl, 2); %figure4
                y_normctrl_std(:,i) = std(yval_normctrl,[],2);

            end %loop timepoints

            for s = 1:2
                maxrange(:,s) = max(yvaluesnormtoTP0(:,s)) - min(yvaluesnormtoTP0(:,s));  %maximal range of final responses
                meanresp(:,s) = mean(yvaluesnormtoTP0(:,s),1); %mean of final responses
                stdevresp(:,s) = std(yvaluesnormtoTP0(:,s),[],1); %stdev of final responses
                coeffvarresp(:,s) = stdevresp(:,s)/meanresp(:,s); %coefficient of variation of final responses
                aucresp(:,s) = trapz(timepoints, yvaluesnormtoTP0(:,s)); %AUC finalresponse
                aucdivby24(:,s) = aucresp(:,s)/24; %AUC/24
                aocresp(:,s) = trapz(timepoints, (1-yvaluesnormtoTP0(:,s))); %AOC
                aucaocratio(:,s) = aocresp(:,s)/aucresp(:,s); %Ratio AOC/AUC

                ctrl_maxrange(:,s) = max(ctrlvaluesnormtoTP0(:,s)) - min(ctrlvaluesnormtoTP0(:,s));
                ctrl_meanresp(:,s) = mean(ctrlvaluesnormtoTP0(:,s),1);
                ctrl_stdevresp(:,s) = std(ctrlvaluesnormtoTP0(:,s),[],1);
                ctrl_coeffvarresp(:,s) = ctrl_stdevresp(:,s)/ctrl_meanresp(:,s);
                ctrl_aucresp(:,s) = trapz(timepoints, ctrlvaluesnormtoTP0(:,s));
                ctrl_aucdivby24(:,s) = ctrl_aucresp(:,s)/24;
                ctrl_aocresp(:,s) = trapz(timepoints, (1-ctrlvaluesnormtoTP0(:,s)));
                ctrl_aucaocratio(:,s) = ctrl_aocresp(:,s)/ctrl_aucresp(:,s);

                meaninh(s,d) = mean(inhibition(:,s),1);
                stdinh(s,d) = std(inhibition(:,s),[],1);

            end

            ydata_finalTP(:,d) = mean(yval_finalTP,2);
            ydata_finalTP_std(:,d) = std(yval_finalTP,[],2);

            ydata_normtoTP0(:,d) = mean(yvaluesnormtoTP0,2); %figure5.1
            ydata_normtoTP0_std(:,d) = std(yvaluesnormtoTP0,[],2);

            meanx = mean(xdata, 2);

            valuestoclear2 = {'xdata';'yval_raw';'yval_smooth';'yval_normx0';'ctrl_smooth';'ctrlval_normx0';
                'yval_norm_ToT';'ctrl_norm_ToT';'yval_normctrl';'yval_finalTP';'ctrlval_finalTP';
                'yvaluesnormtoTP0';'ctrlvaluesnormtoTP0'};

            clear(valuestoclear2{:});

            %% save final values in a single matrix

            valuestosave_treat = {maxrange,meanresp,stdevresp,coeffvarresp,aucresp,aucdivby24,aucaocratio};

            for vv = 1:numel(valuestosave_treat)
                finalvalues(vv,d) = mean(valuestosave_treat{vv},2);
                finalvalues_std(vv,d) = std(valuestosave_treat{vv},[],2);
            end

            if experiment == 20221107 && d == 9 %NaCl control for Cisplatin
                plus = 3;
            elseif experiment == 20221107 && d == 8 %DMSO control for Adavosertib
                plus = 2;
            else
                plus = 1; %DMSO control for all other drugs
            end

            valuestosave_ctrl = {ctrl_maxrange,ctrl_meanresp,ctrl_stdevresp,ctrl_coeffvarresp,ctrl_aucresp,ctrl_aucdivby24,ctrl_aucaocratio};

            for vv = 1:numel(valuestosave_treat)
                finalvalues(vv,dd+plus) = mean(valuestosave_ctrl{vv},2);
                finalvalues_std(vv,dd+plus) = std(valuestosave_ctrl{vv},[],2);
            end

            %% Figure 1: Raw Data
            fig1 = figure(1);

            for j = 1:ii

                x = transpose(meanx);
                y1 = transpose(ydata_raw(:,j));
                error1 = transpose(ydata_raw_std(:,j));
                ctrl1 = transpose(ctrl_raw(:,j));
                ctrlerror1 = transpose(ctrl_raw_std(:,j));

                patch([x fliplr(x)], [(y1-error1)  (fliplr(y1+error1))], color{j}, ...
                    'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold on
                plot(x,y1,'LineWidth',2,'LineStyle','-','Color',color{j}); hold on

                patch([x fliplr(x)], [(ctrl1-ctrlerror1)  (fliplr(ctrl1+ctrlerror1))], color{j}, ...
                    'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold all
                plot(x,ctrl1,'LineWidth',2,'LineStyle','--','Color',color{j},'HandleVisibility','off'); hold on

            end

            % Create labels and title
            ax = gca;
            if a == 1
                ylabeltext1 = append(yaxisnames{a},' [%]');
                yticks(0:20:100); ylim([0 110]);
            else
                ylabeltext1 = append(yaxisnames{a});
                ax.YAxis.Exponent = 3;
            end

            mainTextOpt = {'FontSize',22,'FontName','Helvetica Neue','FontWeight','normal'};
            xlabel('Time (h)',mainTextOpt{:});
            ylabel(ylabeltext1,mainTextOpt{:});
            titletext = append(cellline{c},' ',drug{d});
            title(titletext,mainTextOpt{:});

            % Set the remaining axes and box properties
            grid on;
            xticks(0:24:192);
            if experiment == 20220228 || experiment == 20221107 || experiment == 20211001 || experiment == 20211027
                xticklabels({'-48','-24','0','24','48','72','96','120','144'});
            elseif experiment == 20220215
                xticklabels({'-24','0','24','48','72','96','120'});
            elseif experiment == 20220209
                xticklabels({'0','24','48','72','96','120'});
            end

            axOpt = {'XLimitMethod','padded','YLimitMethod','padded','linewidth',1.5,'YGrid','on','XGrid','off','Box','on','Color','none','FontSize',22};
            set(ax,axOpt{:});

            % add legend
            lgdOpt = {'FontSize',14,'Orientation','horizontal','Location','north', ...
                'FontWeight','normal','EdgeColor','none','Color','#f5f5f5'};
            legendentries = num2str(timepoints);
            lgd = legend(legendentries);
            set(lgd,lgdOpt{:});

            title(lgd,'Time of Day (h)','FontWeight','normal','FontName','Helvetica Neue','FontSize',16);

            %save figure1
            %             filename1 = append(date,'_ToD_',cellline{b},'_',drug{d},'_',channel{a},'_fig1_RawData');
            %             saveas(fig1, [filename1 '.svg']);

            clf(fig1,'reset');

            %% Figure 2: Data normalized to x=0
            fig2=figure(2);

            for k = 1:ii

                y2 = transpose(ydata_normx0(:,k));
                error2 = transpose(ydata_normx0_std(:,k));
                ctrl2 = transpose(ctrl_normx0(:,k));
                ctrlerror2 = transpose(ctrl_normx0_std(:,k));

                patch([x fliplr(x)], [(y2-error2)  (fliplr(y2+error2))], color{k}, ...
                    'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold all
                plot(x,y2,'LineWidth',2,'LineStyle','-','Color',color{k}); hold on

                patch([x fliplr(x)], [(ctrl2-ctrlerror2)  (fliplr(ctrl2+ctrlerror2))], color{k}, ...
                    'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold all
                plot(x,ctrl2,'LineWidth',2,'LineStyle','--','Color',color{k},'HandleVisibility','off'); hold on

            end

            %Create labels and title
            ylabeltext2 = append('Normalized ',yaxisnames{a});
            xlabel('Time (h)',mainTextOpt{:});
            ylabel(ylabeltext2,mainTextOpt{:});
            title(titletext,mainTextOpt{:});

            ax = gca;
            grid on;

            xticks(0:24:192);
            if experiment == 20220228 || experiment == 20221107 || experiment == 20211001 || experiment == 20211027
                xticklabels({'-48','-24','0','24','48','72','96','120','144'});
            elseif experiment == 20220215
                xticklabels({'-24','0','24','48','72','96','120'});
            elseif experiment == 20220209
                xticklabels({'0','24','48','72','96','120'});
            end
            %             ylim([0.5 Inf]);
            set(ax,axOpt{:});

            %add legend
            lgd = legend(legendentries);
            set(lgd,lgdOpt{:});
            title(lgd,'Time of Day (h)','FontWeight','normal','FontName','Helvetica Neue','FontSize',16);

            xline(treat1, 'color', [1 0 0], 'LineWidth',10,'LineStyle','-','HandleVisibility','off','Alpha',0.1);
            xline(treat2, 'color', [1 0 0], 'LineWidth',10,'LineStyle','-','HandleVisibility','off','Alpha',0.1);

            %save figure 2
            %             filename2 = append(date,'_ToD_',cellline{b},'_',drug{d},'_',channel{a},'_fig2_','Normx0');
            %             saveas(fig2, [destination filename2 '.svg']);

            clf(fig2,'reset');

            %% Figure 3: Data normalized to time of treatment

            fig3 = figure(3);

            for h = 1:ii

                y3 = transpose(ydata_normToT(:,h));
                error3 = transpose(ydata_normToT_std(:,h));
                ctrl3 = transpose(ctrl_normToT(:,h));
                ctrlerror3 = transpose(ctrl_normToT_std(:,h));

                length = size(y3,2);
                x3 = (1:intervals:(length*intervals));

                patch([x3 fliplr(x3)], [(y3-error3)  (fliplr(y3+error3))], color{h}, ...
                    'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold all
                plot(x3,y3,'LineWidth',2,'LineStyle','-','Color',color{h});

                patch([x3 fliplr(x3)], [(ctrl3-ctrlerror3)  (fliplr(ctrl3+ctrlerror3))], color{h}, ...
                    'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold all
                plot(x3,ctrl3,'LineWidth',2,'LineStyle','--','Color',color{h},'HandleVisibility','off'); hold on

            end

            %Create labels and title
            xlabel('Time (h)',mainTextOpt{:});
            ylabel(ylabeltext2,mainTextOpt{:});
            title(titletext,mainTextOpt{:});

            ax = gca;
            grid on;

            xticks(0:24:192);

            %             ylim([0.5 Inf]);
            set(ax,axOpt{:});

            %add legend
            lgd = legend(legendentries);
            set(lgd,lgdOpt{:});
            title(lgd,'Time of Day (h)','FontWeight','normal','FontName','Helvetica Neue','FontSize',16);

            %save figure 3
            %             filename3 = append(date,'_ToD_',cellline{b},'_',drug{d},'_',channel{a},'_fig9_','NormToT_noctrl');
            %             saveas(fig3, [destination filename3 '.svg']);

            clf(fig3,'reset');

            %% Figure 4: Relative growth to respective control
            fig4=figure(4);

            for l = 1:ii

                y4 = transpose(y_normctrl(:,l));
                error4 = transpose(y_normctrl_std(:,l));

                length = size(y4,2);
                x4 = (1:intervals:(length*intervals));

                patch([x4 fliplr(x4)], [(y4-error4)  (fliplr(y4+error4))], color{l}, ...
                    'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold all
                plot(x4,y4,'LineWidth',2,'LineStyle','-','Color',color{l}); hold on

            end

            %Create labels and title
            ylabeltext4 = append('Relative ',yaxisnames{a});
            xlabel('Time (h)',mainTextOpt{:});
            ylabel(ylabeltext4,mainTextOpt{:});
            title(titletext,mainTextOpt{:});

            ax = gca;
            grid on;
            xticks(0:24:192);
            set(ax,axOpt{:});

            %add legend
            lgd = legend(legendentries);
            set(lgd,lgdOpt{:});
            title(lgd,'Time of Day (h)','FontWeight','normal','FontName','Helvetica Neue','FontSize',16);

            %insert line at y=1 (=control)
            yline(1, 'color', [0.4, 0.4, 0.4], 'LineWidth',2,'LineStyle','--','HandleVisibility','off');

            %save figure4
            %             filename4 = append(date,'_ToD_',cellline{b},'_',drug{d},'_',channel{a},'_fig4_NormCtrl');
            %             saveas(fig4, [destination filename4 '.svg']);

            clf(fig4,'reset');

            %% Figure 5: ToD-Response Curves

            if d == 1
                rrcc = 2;
            elseif experiment == 20221107 && d > 7
                rrcc = 2;
            else
                rrcc = 1;
            end

            for rc = 1:rrcc

                x_RC = timepoints;

                if rc == 1 %treated
                    y5 = ydata_normtoTP0(:,d);
                    error5 = ydata_normtoTP0_std(:,d);
                    fignr = 'fig5';
                    drugname = drug{d};
                    method = '';
                    ylineval = 1;
                elseif rc == 2 %control
                    y5 = ctrl_normtoTP0(:,d);
                    error5 = ctrl_normtoTP0_std(:,d);
                    fignr = 'fig6';
                    if d == 1
                        drugname = 'DMSO';
                    elseif d == 8
                        drugname = 'DMSO-Ada';
                    elseif d == 9
                        drugname = 'NaCl';
                    end
                    method = '';
                    ylineval = 1;
                end

                titletext_RC = append(cellline{c},' ',drugname,' ',method);

                fig5 = figure(5);

                errorbar(x_RC,y5,error5,'-o','MarkerSize',12,'MarkerFaceColor',[0 0 0],'LineWidth',2,'Color',[0 0 0]);

                xlabel('Time of Day (h)',mainTextOpt{:});
                ylabeltext5 = append(yaxisnames{a},' relative to ToD-0');
                ylabel(ylabeltext5,mainTextOpt{:});
                title(titletext_RC,mainTextOpt{:});

                %insert line at y=1 (=control)
                yline(ylineval, 'color', [0.4, 0.4, 0.4], 'LineWidth',2,'LineStyle','--');

                ax = gca;
                grid on;
                xlim([-1 25]);
                if experiment == 20211001 || experiment == 20211027
                    xticks(0:3:24);
                    xticklabels({'0','3','6','9','12','15','18','21','24'});
                else
                    xticks(0:4:24);
                    xticklabels({'0','4','8','12','16','20','24'});
                end

                %yticks(0:0.5:2); ylim([0 2]);
                set(ax,axOpt{:});

                %save figure5
                %                 filename5 = append(date,'_ToD_',cellline{b},'_',drugname,'_',channel{a},'_',fignr,'_RC_4d',method);
                %                 saveas(fig5, [destination filename5 '.svg']);

                clf(fig5,'reset');
            end

            %% Figure 6: ToD-Response Curve, spline smoothed + unnormalized to ToD-0

            if d == 1
                rrcc = 2;
            elseif experiment == 20221107 && d > 7
                rrcc = 2;
            else
                rrcc = 1;
            end

            for rc = 1:rrcc

                x_RC = timepoints;

                if rc == 1
                    y6 = response_finalTP(:,d);
                    error6 = response_finalTP_std(:,d);
                    fignr = 'fig7';
                    drugname = drug{d};
                    method = 'spline-unnormalized';
                elseif rc == 2
                    y6 = ctrl_finalTP(:,d);
                    error6 = ctrl_finalTP_std(:,d);
                    fignr = 'fig8';
                    if d == 1
                        drugname = 'DMSO';
                    elseif d == 8
                        drugname = 'DMSO-Ada';
                    elseif d == 9
                        drugname = 'NaCl';
                    end
                    method = 'spline-unnormalized';
                end

                titletext_RC = append(cellline{c},' ',drugname,' ',method);

                fig6 = figure(6);

                splinefit = fit(timepoints,y6,'smoothingspline','SmoothingParam',0.7);

                errorbar(x_RC,y6,error6,'o','MarkerSize',12,'MarkerFaceColor','k','LineWidth',2,'Color','k'); hold on

                % Evaluate spline fit within the range of timepoints
                xRange = linspace(min(timepoints), max(timepoints), 1000); % Adjust the number of points as desired
                yfitinterpolate = feval(splinefit, xRange);
                h1 = plot(xRange, yfitinterpolate);
                hold off

                set(h1,'LineWidth',2.5,'color','k','HandleVisibility','off');

                xlabel('Time of Day (h)',mainTextOpt{:});
                ylabeltext5 = append('Normalized ',yaxisnames{a},' Day 4');
                ylabel(ylabeltext5,mainTextOpt{:});
                title(titletext_RC,mainTextOpt{:});

                ax = gca;
                grid on;
                xlim([-1 25]);
                if experiment == 20211001 || experiment == 20211027
                    xticks(0:3:24);
                    xticklabels({'0','3','6','9','12','15','18','21','24'});
                else
                    xticks(0:4:24);
                    xticklabels({'0','4','8','12','16','20','24'});
                end

                %yticks(0:0.5:2); ylim([0 2]);
                set(ax,axOpt{:});

                %save figure6
                %                     filename6 = append(date,'_ToD_',cellline{b},'_',drugname,'_',channel{a},'_',fignr,'_RC_4d_',method);
                %                     saveas(fig6, [destination filename6 '.svg']);

                clf(fig6,'reset');


                %% Figure 7: ToD-Response Curve, spline smoothed + normalized to ToD-0

                if d == 1
                    rrcc = 2;
                elseif experiment == 20221107 && d > 7
                    rrcc = 2;
                else
                    rrcc = 1;
                end

                for rc = 1:rrcc

                    x_RC = timepoints;

                    if rc == 1
                        y7 = ydata_normtoTP0(:,d);
                        error7 = ydata_normtoTP0_std(:,d);
                        fignr = 'fig9';
                        drugname = drug{d};
                        method = 'spline-normalized';
                    elseif rc == 2
                        y7 = ctrl_normtoTP0(:,d);
                        error7 = ctrl_normtoTP0_std(:,d);
                        fignr = 'fig10';
                        if d == 1
                            drugname = 'DMSO';
                        elseif d == 8
                            drugname = 'DMSO-Ada';
                        elseif d == 9
                            drugname = 'NaCl';
                        end
                        method = 'spline-normalized';
                    end

                    titletext_RC = append(cellline{c},' ',drugname,' ',method);

                    fig7 = figure(7);

                    splinefit = fit(timepoints,y7,'smoothingspline','SmoothingParam',0.7);

                    errorbar(x_RC,y7,error7,'o','MarkerSize',12,'MarkerFaceColor','k','LineWidth',2,'Color','k'); hold on

                    % Evaluate spline fit within the range of timepoints
                    xRange = linspace(min(timepoints), max(timepoints), 1000); % Adjust the number of points as desired
                    yfitinterpolate = feval(splinefit, xRange);
                    h1 = plot(xRange, yfitinterpolate);
                    hold off

                    set(h1,'LineWidth',2.5,'color','k','HandleVisibility','off');

                    xlabel('Time of Day (h)',mainTextOpt{:});
                    ylabeltext5 = append('Normalized ',yaxisnames{a},' Day 4');
                    ylabel(ylabeltext5,mainTextOpt{:});
                    title(titletext_RC,mainTextOpt{:});

                    ax = gca;
                    grid on;
                    xlim([-1 25]);
                    if experiment == 20211001 || experiment == 20211027
                        xticks(0:3:24);
                        xticklabels({'0','3','6','9','12','15','18','21','24'});
                    else
                        xticks(0:4:24);
                        xticklabels({'0','4','8','12','16','20','24'});
                    end

                    %yticks(0:0.5:2); ylim([0 2]);
                    set(ax,axOpt{:});

                    %save figure7
                    %                     filename7 = append(date,'_ToD_',cellline{b},'_',drugname,'_',channel{a},'_',fignr,'_RC_4d_',method);
                    %                     saveas(fig7, [destination filename8 '.svg']);

                    clf(fig7,'reset');

                end

            end
        end %Loop d drug

        outputsheets = {'finalvaluesmean';'finalvaluesstd';'responses';'stdresponses'};
        valuename = {'maxmimum range';'mean';'stdev';'coeffvar';'AUC';'AUC/24';'AOC/AUC'};

        if experiment ~= 20221107
            ydata_normtoTP0(:,end+1) = ctrl_normtoTP0(:,1);
            ydata_normtoTP0_std (:,end+1) = ctrl_normtoTP0_std(:,1);
        else
            gg = [1,8,9]; %add responses of 'DMSO','DMSO-Ada','NaCl' to array
            for g = 1:3
                ydata_normtoTP0(:,dd+g) = ctrl_normtoTP0(:,gg(:,g));
                ydata_normtoTP0_std(:,dd+g) = ctrl_normtoTP0_std(:,gg(:,g));
            end
        end

        input = {finalvalues;finalvalues_std;ydata_normtoTP0;ydata_normtoTP0_std};
        outputfile = append(date,'_ToD_Results_4d_',cellline{c},'.xlsx');

        for s = 1:4
            
            if s == 1 || s == 2
                t_firstcol = cell2table(valuename);
            else
                t_firstcol = array2table(timepoints);
            end

            finalvaluestosave = input{s};

            if experiment == 20220228 && c == 3 && s < 3
                finalvaluestosave(:,3) = NaN;
            end

            if experiment == 20221107
                varname = [(append(drug));'DMSO';'DMSO-Ada';'NaCl'];
            elseif experiment == 20211027 || experiment == 20211001
                varname = [(append(drug));'NaCl'];
            else
                varname = [(append(drug));'DMSO'];
            end

            outputsheet = append(outputsheets{s},'_',channel{a});
            t_finalvaluestosave = array2table(finalvaluestosave,"VariableNames",varname);
            table = [t_firstcol,t_finalvaluestosave];
            writetable(table,outputfile,"sheet",outputsheet);

            clear t_firstcol
            clear finalvaluestosave
            clear table

        end

        valuestoclear3 = {'meanx';'ydata_raw';'ydata_raw_std';'ctrl_raw';'ctrl_raw_std';'ydata_normx0';'ydata_normx0_std';
            'ctrl_normx0';'ctrl_normx0_std';'y_normctrl';'y_normctrl_std';'ydata_normtoTP0';'ydata_normtoTP0_std';'ctrl_normtoTP0';'ctrl_normtoTP0_std'};
        clear(valuestoclear3{:});

    end %loop c cell lines
end  %loop a channel
end %function
