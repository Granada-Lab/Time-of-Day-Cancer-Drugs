function polarhistogram_ToD_max_min_benefit_min_diff_cellline(file,excelsheets_times,channel,celllines)

%Carolin Ector, 23.08.2023

%Function creates polarhistograms/circular histograms to evaluate the times of:
% - maximum benefit (highest foldchange between cancer cell line and mcf10a with higher toxicity against cancer cell line as opposed to mcf10a)
% - minimum benefit (highest foldchange between cancer cell line and mcf10a with higher toxicity against mcf10a as opposed to cancer cell line)
% - minimum difference (lowest foldchange/difference in responses between cancer cell line and mcf10a)

%Time-of-Day-Cancer-Drugs Manuscript Fig. 4i

%input: stored in "ToD_analysis_workspace.mat"
% file: excel file where maximum and minumum benefit foldchanges and times are stored
% excelsheets: sheet names for the time of the day (h)
% celllines: names of the cell lines being analysed
% channel: either brightfield (confluency) or fluorescent (cell number) channel from live-imaging, in the manuscript we focus on the fluorescent channel
% drugs: names of drugs the cell lines have been treated with

%define variables
titles = {'ToD Max and Min Benefit';'ToD Minimum Difference'};
colorall1 = {[0.12 0.82 0.37];[1 0 0];'none'};
colorall2 = {'none';'none';[0 0 0]};
edgecolor1 = {[0 0 0];[0 0 0];'none'};
edgecolor2 = {'none';'none';[0 0 0]};

legendentries = {'maximum benefit';'minimum benefit'};

cc = numel(celllines);
ii = numel(excelsheets_times);

for a = 1:2 %loop a channel

    for c = 1:cc %loop c cellline

        fig = figure;
        fig.Position = [787,627,1401,533];
        sgtitle(celllines{c},'FontSize',14,'FontName','Helvetica Neue','FontWeight','bold');

        for i = 1:ii %loop input

            %load data
            inputsheet = append(excelsheets_times{i},channel{a});
            [alldata] = readmatrix(file,'sheet',inputsheet);
            alldata(:,[1,2,9]) = [];
            %save data from all drugs per cell line
            data(:,i) = alldata(:,c);

        end

        % Find indices where the value is 24 
        indices = data == 24;
        % Replace those values with 23.5 --> so that they fall into 22-24 hour bin and not in 24-2 hour bin in the polarhistogram
        data(indices) = 23.5;

        for k = 1:2 %loop k: create two polarhistograms side-by-side

            if k == 1
                subplot(1,2,1) %overlay max and min benefit times with different colors
            else
                subplot(1,2,2) %overlay max and min benefit times + min difference, only show min difference
            end

            for j = 1:3

                %calculate number of samples the polarhistogram data is based on
                n = size(rmmissing(data(:,j),1));
                titletext = append(titles{k},' n=',num2str(n(:,1)));

                if k == 1
                    color = colorall1{j};
                    edgecolor = edgecolor1{j};
                else
                    color = colorall2{j};
                    edgecolor = edgecolor2{j};
                end

                %plot data from all drugs per cell line
                alpha = fliplr(deg2rad(rmmissing(data(:,j)/ 24 * 360)));
                polarhistogram(alpha,'Normalization','probability','BinEdges',deg2rad(0:30:360),'FaceColor',color,'EdgeColor',edgecolor,'FaceAlpha',0.4); hold all;

                %create legend and fine-tune plot
                if k == 1
                    legend(legendentries,'Location','northeastoutside');
                end

                ax = gca;
                ax.PositionConstraint = "innerposition";
                title(titletext,'FontSize',22,'FontName','Helvetica Neue','FontWeight','normal');
                set(ax,'linewidth',3,'Color','none','FontSize',22,'FontName','Helvetica Neue');

                % Set custom theta tick labels
                thetaticks([0:30:360]);
                thetaticklabels({'0/24','2','4','6','8','10','12','14','16','18','20','22','24'});
                set(gca, 'ThetaZeroLocation', 'top');
                set(gca, 'RAxisLocation', 330);

            end %loop input
        end %loop polarhistograms

        %save figure
        filename = append('ToD_plots/',channel{a},'/polarhistogram_max_min_times/polarhistogram_',celllines{c},'_',channel{a});
        saveas(fig, [filename '.svg']);

    end %loop c cellline

end %loop a channel

end %function