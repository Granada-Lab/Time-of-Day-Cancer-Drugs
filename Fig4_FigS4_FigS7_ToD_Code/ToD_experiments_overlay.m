function ToD_experiments_overlay(c,experiments,drugname,cellline)

%Carolin Ector
%based on 'ToD_U2OS_overlay.m' (20.03.2024)
%modified on 25.03.2024 to overlay different MCF7 experiments

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
channel = {'Conf';'CellNr'}; %confluency = brightfield channel, cell number = red fluorescent channel from live-imaging
yaxisnames = {'Confluency (%)';'Cell Number'};

colorarray = lines(numel(experiments));

%% Load and normalize data

for a = 2:numel(channel) %loop a channel

    sheet1 = append('responses_',channel{a});
    sheet2 = append('stdresponses_',channel{a});

    for d = 1:numel(drugname) %loop d drug

        figure('Position',[1000,918,672,420]);

        for e = 1:numel(experiments) %loop e experiments

            if c == 1 %MCF7
                if d <= 4
                    exp = append(experiments{e},'1');
                    col = d;
                elseif d >= 5
                    exp = append(experiments{e},'2');
                    col = d-4;
                end
            else
                exp = append(experiments{e});
                col = d;
            end

            % Load excel file where data is stored
            inputfile = append('ToD_results/',exp,'_ToD_Results_4d_',cellline,'.xlsx');

            [data] = readmatrix(inputfile{1},'Sheet',sheet1);
            [stdev] = readmatrix(inputfile{1},'Sheet',sheet2);

            timepoints = data(:,1);
            responses = data(:,col+1);
            err = stdev(:,col+1);
            coeffvar(e,d) = mean(err,'all')/mean(responses,"all");

            %% spline smoothing

            %Spline interpolation provides a smooth curve that passes through the given data points.
            splinefit = fit(timepoints,responses,'smoothingspline','SmoothingParam',0.7);

            % Evaluate spline fit within the range of timepoints
            xRange = linspace(min(timepoints), max(timepoints), 1000); % Adjust the number of points as desired
            yfitinterpolate(:,e) = feval(splinefit, xRange);
            ToDMR_spline(e,d) = max(yfitinterpolate(:,e))-min(yfitinterpolate(:,e)); %save ToDMR value from spline smooth

            h1 = plot(xRange, yfitinterpolate(:,e)); hold all
            % fine-tune plots
            set(h1,'LineWidth',2,'color',colorarray(e,:),'HandleVisibility','off');
            errorbar(timepoints,responses,err,'o','MarkerFaceColor',colorarray(e,:),'Color',colorarray(e,:));

            clear yfitinterpolate

        end %loop e experiments

        % insert line at y=1 (=control)
        yline(1, 'color', [0.4, 0.4, 0.4], 'LineWidth',2,'LineStyle','--','HandleVisibility','off');

        % add legend
        hLeg = legend(experiments);
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
        filename = append('ToD_plots/overlay_MCF7_experiments/ToD_',drugname{d},'_',channel{a},'_splinesmooth');
        saveas(h1, [filename '.svg']);

    end %loop d drug

    %close all

    %     for d2 = 1:numel(drugname) %loop d2 drug
    %
    %         fig = figure('Visible','off'); %plot raw ToDMR values per drug
    %
    %         x = ToDMR_spline(:,d2);
    %         figurtext = 'raw_ToDMR_spline_raw';
    %         ylabeltext = 'ToD_{MR}';
    %
    %         bar(x,'BarWidth',0.7,'FaceColor',[0.9,0.9,0.9],'LineWidth',0.9);
    %         xticklabels(labels);
    %         title(drugname{d2});
    %         ylabel(ylabeltext);
    %         ax = gca;
    %         axOpt = {'linewidth',0.8,'YGrid','off','XGrid','off','Box','on','Color','none','FontSize',14};
    %         set(ax,axOpt{:});
    %
    %         xaxisproperties= get(gca, 'XAxis');
    %         xaxisproperties.TickLabelInterpreter = 'none';
    %
    %         %save figure
    %         filename2 = append('ToD_plots/overlay_MCF7_experiments/ToD_',drugname{d},'_',channel{a},'_splinesmooth_bar_',figurtext);
    %         saveas(fig, [filename2 '.svg']);
    %
    %         clear x
    %         clear figuretext
    %         clear ylabeltext
    %         clear labels
    %
    %     end  %loop d2 drug


    fig2 = figure; %plot raw ToDMR values per drug

    heatmapylabels = drugname;

    Y = ToDMR_spline;
    titletext = 'ToD_{MR}';

    h1 = heatmap(experiments,heatmapylabels,Y');
    h1.Title = titletext;
    h1.FontSize = 14;
    h1.CellLabelFormat = '%0.1g';
    %         s = struct(h1);
    %         s.XAxis.TickLabelRotation = 45;  %horizontal orientation of the bar plot
    colormap("parula");

    clear Y
    clear titletext

    %save figure
    filename3 = append('ToD_',channel{a},'_splinesmooth_heatmap');
    saveas(fig2, [filename3 '.svg']);

    meancoeffvarall = mean(coeffvar,2);
    disp(meancoeffvarall)
    disp(channel{a})
    
%     %save ToDMR values (absolute and relative to U2OS WT) from spline smooth
%     input = {ToDMR_spline};
%     outputsheetnames = {'ToDMR_spline_'};
%     rownames = cell2table(experiments);
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