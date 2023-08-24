% Normalize Excel Data and plot ToD Response Curves

function expfit_growthrates_ToD_manuscript_selected_celllines(cellline,colorall)

%Fits an exponential function to growth curves, plate by plate
%INPUT
% cellline = all cellline as string in cell
%OUTPUT
% k-value from fit
% Ymax value from growth curve
% Relative Growth Factor(RGF) = (relative k) * (relative ymax)

%% Define variables
cc = numel(cellline);
metric = {'CellNr';'Conf'};
yaxisnames = {'Cell Number';'Confluence'};
initialvaluetext = {'Cell Number_{t0} = ';'Confluence_{t0} = '};

% Load excel file where data is stored
inputfile = 'Growth_Data_48_well_Selection_Manuscript_v2';

%destination output (figures, ecxel sheets)
destination = append('selected_celllines_manuscript/');

for m = 1:1 %loop metric

    %create excel file where data will be stored

    outputfile = append(destination,'Results_Growth_Data_48_well_Selection_Manuscript_',metric{m},'_v2.xlsx');
    yaxisname = append('Normalized ',yaxisnames{m});

    sheet1 = append('Data_',metric{m});
    sheet2 = append('Error_',metric{m});

    [data] = xlsread(inputfile,sheet1);
    [error] = xlsread(inputfile,sheet2);

    fig1 = figure;
    fig1.Position = [1,1,2560,1361];

    fitresult = cell(cc,1); %this is saved to the workspace in an answer array and it contains the equation and the calculated values
    gof = struct( 'sse', cell(cc,1), 'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] ); %gof is somehow not saved to the workspace

    for c = 1:cc %loop a cell lines

        if ismember(c, [1, 7, 8, 11, 12])
            time = (0:1.5:96)';
        else
            time = (0:2:96)';
        end

        %take data from the first 4 days to make different experiments with different recording lengths comparable
        endrow=size(time,1);

        % load, smooth and normalize data
        rawdata = data(1:endrow,c);
        smootheddata = smooth(rawdata,0.35,'rloess');
        normalizeddata = smootheddata./smootheddata(1,:);

        stdev = error(1:endrow,c);
        smoothedstdev = smooth(stdev,0.35,'rloess');
        normalizedstdev = smoothedstdev./smootheddata(1,:);

        normalizeddata_all{m,c} = normalizeddata;
        stdnormalizeddata_all{m,c} = normalizedstdev;
        time_all{c} = time;

        smootheddata_all{m,c} = smootheddata;
        stdsmootheddata_all{m,c} = smoothedstdev;

        % calculate doubling times
        yend = median(normalizeddata(end-3:end,:));
        SE_yend = median(normalizedstdev(end-3:end,:));
        ystart = median(normalizeddata(1:3,:));
        SE_ystart = median(normalizedstdev(1:3,:));
        DoublingTime{1,c} = 96*log(2)/log(yend/ystart);
        DoublingTime{2,c} = (96 * log(2) / (log(yend / ystart))^2) * sqrt((SE_yend / yend)^2 + (SE_ystart / ystart)^2);

        % save initial cell count
        initialvalue{c} = median(smootheddata(1:3,:));

        %% Fit in a loop

        [xData, yData] = prepareCurveData(time, normalizeddata);

        %exponential curve fit
        fitfunc = 'expfit';
        equation = 'a*exp(x*k)';

        %         fitfunc = 'logfit';
        %         equation = 'Ymax/(1+exp(-k*(x-x0)))';

        % Set up fittype and options.
        ft = fittype( equation, 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Lower = [0.9 -0.1]; % [Ymin k]
        opts.Robust = 'on';
        opts.StartPoint = [1 0.5]; % [Ymin k]
        opts.Upper = [1.1 0.1]; % [Ymin k]

        % Fit model to data.
        [fitresult{c}, gof(c)] = fit( xData, yData, ft, opts );
        vectcoeff=coeffvalues(fitresult{c});
        % assign k values to allkvalues variable which is saved to workspace
        Rsq{c} = (gof(c,1).rsquare);
        growthrates{1,c} = vectcoeff(2);
        confBounds = confint(fitresult{c});
        growthrates{2,c} = confBounds(1,2); %lower confidence bounds growthrate
        growthrates{3,c} = confBounds(2,2); %upper confidence bounds growthrate

        if mod(cc, 2) == 0
            nr = cc/2;
        else
            nr = (ceil(cc/2)*2)/2;
        end

                h=subplot(2,nr,c);
                hold all
                x = transpose(xData);
                y = transpose(yData);
                e = transpose(normalizedstdev);
                patch([x, fliplr(x)], [y - e, fliplr(y + e)], 'b','FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'HandleVisibility', 'off');
                plot( fitresult{c}, xData, yData );
        
                if m == 1
                    parameters=strcat(['Cell Number_{t0} = ',num2str(initialvalue{c},'%6.1f'),newline,'\color{black}Doubling Time = ',num2str(DoublingTime{c},'%6.1f'),' hours',newline,'\color{red}k = ',num2str(vectcoeff(2),'%6.3f'),newline,'\color{black}R^2 = ',num2str(gof(c,1).rsquare,'%6.3f')]);           
                else
                    parameters=strcat(['Confluence_{t0} = ',num2str(initialvalue{c},'%6.1f'),' %',newline,'\color{black}Doubling Time = ',num2str(DoublingTime{c},'%6.1f'),' hours',newline,'\color{red}k = ',num2str(vectcoeff(2),'%6.3f'),newline,'\color{black}R^2 = ',num2str(gof(c,1).rsquare,'%6.3f')]);
                end
               
                ylabel(h,[]);
                xlabel(h,[]);
        
                %create legend
                firstlegendentry = cellline{c};
                legend( h, firstlegendentry, parameters, 'Location', 'northoutside', 'Interpreter', 'tex','color','w');
        
                % Set the remaining axes and box properties
                ax = gca;
                grid on;
                xticks(0:24:120);
                set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
                    'XGrid','off','Box','on','Color','none','FontSize',14);

        varstoclear = {'xData';'yData';'rawdata';'smootheddata';'normalizeddata';'stdev';'smoothedstdev';'normalizedstdev';'vectcoeff';'time'};
        clear(varstoclear{:});

    end %loop c cell lines

    %     % Give common xlabel, ylabel and title to your figure
    %     han=axes(fig1,'visible','off');
    %     titletext = append('fitfunc ','[ f(x) = ',equation,' ]',newline,' ');
    %     han.Title.Visible='on'; han.XLabel.Visible='on'; han.YLabel.Visible='on';
    %     xlabel(han,'Time [hours]','FontWeight','bold','FontSize',16);
    %     ylabel(han,yaxisname,'FontWeight','bold','FontSize',16);
    %     title(han,titletext,'FontSize',16);
    %
    %     hold off
    %
    %     %save figure
    %     filetext = append(destination,fitfunc,'_selected_celllines_manuscript_',metric{m},'_v2');
    % %     saveas(fig1, [ filetext, '.svg']);

    %save values
    varnames = cellline;
    valuesgrowthrate = {'growthrate';'confidencebound_low';'confidencebound_high'};
    valuesdoublingtime = {'doublingtime';'standarderror'};
    T_valuesgrowthrate = cell2table(valuesgrowthrate);
    T_valuesdoublingtime = cell2table(valuesdoublingtime);

    T_growthrates = cell2table(growthrates,"VariableNames",varnames);
    T_growthrates2 = [T_valuesgrowthrate,T_growthrates];
    T_Rsq = cell2table(Rsq,"VariableNames",varnames);
    T_DoublingTime = cell2table(DoublingTime,"VariableNames",varnames);
    T_DoublingTime2 = [T_valuesdoublingtime,T_DoublingTime];
    T_initialcellcount = cell2table(initialvalue,"VariableNames",varnames);

    writetable(T_growthrates2,outputfile,'sheet','growthrate');
    writetable(T_Rsq,outputfile,'sheet','Rsq');
    writetable(T_DoublingTime2,outputfile,'sheet','doublingtime');
    writetable(T_initialcellcount,outputfile,'sheet','initialvalue');

end % loop m metric

% for n = 1:2
% 
%     fig2 = figure;
% 
%     yaxisname2 = append('Normalized ',yaxisnames{n});
% 
%     cellline2 = {'MCF10A';'HCC1143';'HCC1937';'HCC38';'MDAMB468';'HCC1806';'SUM149PT';'CAL51';'MDAMB231';'MDAMB436'};
%     dd = numel(cellline2);
%     colorder = [15,1,2,3,4,5,7,9,10,11];
% 
%     for d = 1:dd
% 
%         column = colorder(:,d);
% 
%         x = time_all{1,column};
%         inputdata = normalizeddata_all{n,column};
%         inputstd = stdnormalizeddata_all{n,column};
% 
%         x2 = transpose(x);
%         y2 = transpose(inputdata);
%         error2 = transpose(inputstd);
% 
%         color = str2num(cell2mat(colorall(column,:)));
% 
%         hold all;
%         patch([x2, fliplr(x2)], [y2 - error2, fliplr(y2 + error2)],color,'EdgeColor', 'none','FaceAlpha', 0.1, 'HandleVisibility', 'off');
%         plot(x2,y2,'Color',color,'LineWidth',3);
% 
%     end
% 
%     ylabel(yaxisname2);
%     xlabel('Time [hours]');
% 
%     %create legend
%     legend(cellline2, 'Location', 'NorthWest', 'Interpreter', 'tex','color','w');
% 
%     % Set the remaining axes and box properties
%     ax = gca;
%     grid on;
%     xticks(0:24:120);
%     set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
%         'XGrid','off','Box','on','Color','w','FontSize',20);
% 
%     %save figure
%     filetext2 = append('norm_',metric{n},'_48_well_overlay_all_celllines');
% %     saveas(fig2, [ filetext2, '.svg']);
% 
%     varstoclear = {'x2';'y2';'error2'};
%     clear(varstoclear{:});
% 
% end %loop n
% 
% for g = 1:cc %loop celllines
% 
%     fig3 = figure;
% 
%     x = transpose(time_all{1,g});
% 
% %     y_cellnr = transpose(normalizeddata_all{1,g});
% %     e_cellnr = transpose(stdnormalizeddata_all{1,g});
% % 
% %     y_conf = transpose(normalizeddata_all{2,g});
% %     e_conf = transpose(stdnormalizeddata_all{2,g});
% 
%     y_cellnr = transpose(smootheddata_all{1,g});
%     e_cellnr = transpose(stdsmootheddata_all{1,g});
% 
%     y_conf = transpose(smootheddata_all{2,g});
%     e_conf = transpose(stdsmootheddata_all{2,g});
% 
%     colororder({'r','b'})
% 
%     hold all;
%     yyaxis left
% %     patch([x, fliplr(x)], [y_cellnr - e_cellnr, fliplr(y_cellnr + e_cellnr)],'r','FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'HandleVisibility', 'off');
%     plot(x,y_cellnr,'Color','r','LineWidth',2);
%     patch([x, fliplr(x)], [y_cellnr - e_cellnr, fliplr(y_cellnr + e_cellnr)],'r','FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'HandleVisibility', 'off');
%     
%     ylabel('Cell Number');
%     yyaxis right
% %     patch([x, fliplr(x)], [y_conf - e_conf, fliplr(y_conf + e_conf)],'b','FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'HandleVisibility', 'off');
%     plot(x,y_conf,'Color','b','LineWidth',2);
%     patch([x, fliplr(x)], [y_conf - e_conf, fliplr(y_conf + e_conf)],'b','FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'HandleVisibility', 'off');
%     
%     ylabel('Confluence');
%     yyaxis right
% 
%     xlabel('Time [hours]');
%     title(cellline{g});
% 
%     %create legend
%     %     legend( b, firstlegendentry, mean_parameters, 'Location', 'NorthWest', 'Interpreter', 'tex','color','w');
%     hold off
% 
%     % Set the remaining axes and box properties
%     ax = gca;
%     grid on;
%     xticks(0:24:120);
%     set(ax,'XLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
%         'XGrid','off','Box','on','Color','w','FontSize',16);
% 
%     %save figure
%     filetext3 = append(destination,'overlay_cellnr_conf_rawsmoothed',cellline{g});
%     saveas(fig3, [ filetext3, '.svg']);
% 
%     varstoclear2 = {'x';'y_cellnr';'e_cellnr';'y_conf';'e_conf'};
%     clear(varstoclear2{:});
%     clf(fig3,'reset');
% end

end %function