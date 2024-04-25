function ToD_U2OS_overlay(date,cellline,drug,timepoints)

%Carolin Ector
%based on 'ToD_mcf10a_overlay_benefit_extraction.m' (23.08.2023)
%modofied on 20.03.2024 to overlay U2OS cell lines

%%Function ...:
%%reads Relative ToD Response data (ToD profile) from excel sheet (calculated with function: ToD_growth_plots_parameter_extraction)
%%overlays ToD profiles of U2OS WT cell lines and clock knockout variants (Cry1-sKO, Cry2-sKO, Cry1/2-dKO)

%Time-of-Day-Cancer-Drugs Manuscript Nature Communications Revision #1

%input: stored in "ToD_U2OS_overlay_workspace.mat"
% date: date of the experiment being analyzed
% cellline: names of the cell lines being analysed
% drug: names of drugs used for time-of-day treatments
% timepoints: timepoints of treatment within the circadian cycle

%Define remaining variables
channel = {'Conf'}; %confluency = brightfield channel, cell number = red fluorescent channel from live-imaging
yaxisnames = {'Confluency (%)'};

cc = numel(cellline);
colorarray = lines(cc);

%% Load and normalize data

for a = 1:numel(channel) %loop a channel

    outputfile = 'ToDMR_splinesmooth_U2OS.xlsx';

    sheet1 = append('responses_',channel{a});
    sheet2 = append('stdresponses_',channel{a});

%     drugname = [(append(drug));'DMSO';'NaCl'];

    drugname = drug;

    for d = 1:numel(drug) %loop d drug

        figure('Position',[1000,918,672,420])%,'Visible','off');

        for c = 1:cc %loop c cell lines

            % Load excel file where data is stored
            inputfile = append('ToD_results/20240311_ToD_Results_4d_',cellline{c},'.xlsx');
            [data] = readmatrix(inputfile,'Sheet',sheet1);
            [stdev] = readmatrix(inputfile,'Sheet',sheet2);

            responses = data(:,d+1);
            err = stdev(:,d+1);

            %% spline smoothing

            %Spline interpolation provides a smooth curve that passes through the given data points.
            splinefit = fit(timepoints,responses,'smoothingspline','SmoothingParam',0.7);

            % Evaluate spline fit within the range of timepoints
            xRange = linspace(min(timepoints), max(timepoints), 1000); % Adjust the number of points as desired
            yfitinterpolate(:,c) = feval(splinefit, xRange);
            ToDMR_spline(c,d) = max(yfitinterpolate(:,c))-min(yfitinterpolate(:,c)); %save ToDMR value from spline smooth

            y_bar(d,c) = ToDMR_spline(c,d);

%             if c ~= cc %overlay U2OS cell lines only (exclude SY5Y)
                h1 = plot(xRange, yfitinterpolate(:,c)); hold all
                % fine-tune plots
                set(h1,'LineWidth',2,'color',colorarray(c,:),'HandleVisibility','off');
                errorbar(timepoints,responses,err,'o','MarkerFaceColor',colorarray(c,:),'Color',colorarray(c,:));
%             end

            clear yfitinterpolate

        end %loop c cell lines

        % insert line at y=1 (=control)
        yline(1, 'color', [0.4, 0.4, 0.4], 'LineWidth',2,'LineStyle','--','HandleVisibility','off');

        % add legend
        hLeg = legend(cellline);
        lgdOpt = {'FontSize',16,'Orientation','vertical','Location','northeastoutside', ...
            'FontWeight','normal','EdgeColor','none','Color','#f5f5f5','Interpreter','none'};
        set(hLeg,lgdOpt{:});

        titletext = drugname{d};

        mainTextOpt = {'FontSize',22,'FontName','Helvetica Neue','FontWeight','normal'};
        axOpt = {'XLimitMethod','padded','YLimitMethod','padded','linewidth',1.5,'YGrid','on','XGrid','off','Box','on','Color','none','FontSize',22};

        xlabel('Time of Day (h)',mainTextOpt{:});
        ylabeltext = append('Relative ',yaxisnames{a});
        ylabel(ylabeltext,mainTextOpt{:});
        title(titletext,mainTextOpt{:});

        ax = gca;
        grid on;

        xticks(0:4:24); xlim([-1 25]);
        set(ax,axOpt{:});

        hold off

        %save figure
%         filename = append('ToD_plots/overlay_WT_KO/',date,'_ToD_',drugname{d},'_',channel{a},'_splinesmooth');
        filename = append('ToD_plots/overlay_WT_KO/',date,'_ToD_',drugname{d},'_',channel{a},'_splinesmooth_dKO_only');
        saveas(h1, [filename '.svg']);

    end %loop d drug

    %close all

    %convert ToDMR = 0 to NaN (no data for Adavosertib and Cisplatin for U2OS Cry2-sKO)
    ToDMR_spline(ToDMR_spline==0) = NaN;

    for d2 = 1:numel(drugname) %loop d2 drug

        for c2 = 1:numel(cellline) %loop c2 cellline
            rel_ToDMR(c2,:) = ToDMR_spline(c2,:)./ToDMR_spline(1,:);%calculate relative ToDMR value to U2OS WT
        end

        for b = 1:2 %loop b barplots

            fig = figure('Visible','off'); %plot raw ToDMR values per drug

            if b == 1
                x = ToDMR_spline(:,d2);
                figurtext = 'raw_ToDMR_spline_raw';
                labels = cellline;
                ylabeltext = 'ToD_{MR}';
            else
                x = rel_ToDMR(:,d2)-1;
                x(1,:) = []; %remove U2OS WT (always 1)
                labels = cellline(2:end,:); %remove U2OS WT from label list
                figurtext = 'relative_ToDMR_spline';
                ylabeltext = 'Relative ToD_{MR} to U2OS-WT';
            end

            bar(x,'BarWidth',0.7,'FaceColor',[0.9,0.9,0.9],'LineWidth',0.9);
            xticklabels(labels);
            title(drugname{d2});
            ylabel(ylabeltext);
            ax = gca;
            axOpt = {'linewidth',0.8,'YGrid','off','XGrid','off','Box','on','Color','none','FontSize',14};
            set(ax,axOpt{:});

            xaxisproperties= get(gca, 'XAxis');
            xaxisproperties.TickLabelInterpreter = 'none';

            %save figure
            filename2 = append('ToD_plots/overlay_WT_KO/',date,'_ToD_',drugname{d},'_',channel{a},'_splinesmooth_bar_',figurtext);
            %saveas(fig, [filename2 '.svg']);

            clear x
            clear figuretext
            clear ylabeltext
            clear labels

        end

    end
% 
%     fig2 = figure; %plot raw ToDMR values per drug
% 
% %     heatmapxlabels = {'U2OS-WT';'U2OS-Cry1-sKO';'U2OS-Cry1/2-dKO';'U2OS-Cry2-sKO'};
%     heatmapxlabels = {'U2OS-WT';'U2OS-Cry1-sKO';'U2OS-Cry1/2-dKO';'U2OS-Cry2-sKO';'SY5Y'};
% 
%     for h = 1:2 %loop h heatmaps
% 
%         heatmapylabels = drug;
% %         heatmapylabels(3,:) = [];
% 
%         if h == 1
% %             Y = ToDMR_spline([1:4],[1:4]); %exclude DMSO and NaCl controls, exclude SY5Y
% %             Y = ToDMR_spline([1:4],[1,2,4]); %exclude DMSO and NaCl controls, exclude SY5Y + Ada
%             Y = ToDMR_spline(:,[1:4]); %exclude DMSO and NaCl controls
%             titletext = 'ToD_{MR}';
%         else
% %             Y = rel_ToDMR([1:4],[1:4]); %exclude DMSO and NaCl controls, exclude SY5Y
% %             Y = rel_ToDMR([1:4],[1,2,4]); %exclude DMSO and NaCl controls, exclude SY5Y + Ada
%             Y = rel_ToDMR(:,[1:4]); %exclude DMSO and NaCl controls
%             titletext = 'ToD_{MR} relative to U2OS-WT';
%         end
% 
%         subplot(1,2,h);
%         h1 = heatmap(heatmapxlabels,heatmapylabels,Y');
%         h1.Title = titletext;
%         h1.FontSize = 14;
%         h1.CellLabelFormat = '%0.1g';
% %         s = struct(h1);
% %         s.XAxis.TickLabelRotation = 45;  %horizontal orientation of the bar plot
%         colormap("parula");
% 
%         clear Y
%         clear titletext
%         clear labels
% 
%     end %loop h heatmaps

    %save figure
%     filename3 = append(date,'_ToD_',channel{a},'_splinesmooth_heatmap');
%     filename3 = append(date,'_ToD_',channel{a},'_splinesmooth_heatmap_noAda');
%     filename3 = append(date,'_ToD_',channel{a},'_splinesmooth_heatmap_withSY5Y');
%     saveas(fig2, [filename3 '.svg']);

    fig3 = figure; 

    y_bar(3,:) = [];
    xlabels = drug;
    xlabels(3,:)= [];
    legendentries = {'U2OS-WT';'U2OS-Cry1/2-dKO'};

    bar(y_bar)
    legend(legendentries);
    xticklabels(xlabels);
    ylabeltext = 'ToD_{MR} value';
    ylabel(ylabeltext);
    xaxisproperties= get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'none';
    ax = gca;
    axOpt = {'linewidth',0.8,'YGrid','off','XGrid','off','Box','on','Color','none','FontSize',14};
    set(ax,axOpt{:});
    
%     %save ToDMR values (absolute and relative to U2OS WT) from spline smooth
%     input = {ToDMR_spline,rel_ToDMR};
%     outputsheetnames = {'ToDMR_spline_';'Rel_ToDMR_'};
%     rownames = cell2table(cellline);
% 
%     for i = 1:numel(input)
% 
%         outputsheet = append(outputsheetnames{i},channel{a});
%         table = array2table(input{i},"VariableNames",drugname);
%         finaltable = [rownames,table];
%         writetable(finaltable,outputfile,'sheet',outputsheet);
% 
%         clear table
%         clear outputsheet
%     end

end %loop a channel

end %function
