function polarhistogram_ToD_max_min_benefit_min_diff_drug(file,excelsheets_times,channel,celllines,drugs)

%Carolin Ector, 23.08.2023
%Time-of-Day-Cancer-Drugs Manuscript Fig. 4i,j

%create polarhistograms/circular histograms to evaluate the times of:
% - maximum benefit (highest foldchange between cancer cell line and mcf10a with higher toxicity against cancer cell line as opposed to mcf10a)
% - minimum benefit (highest foldchange between cancer cell line and mcf10a with higher toxicity against mcf10a as opposed to cancer cell line)
% - minimum difference (lowest foldchange/difference in responses between cancer cell line and mcf10a)

%input: stored in ToD_analysis_workspace.mat
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
dd = numel(drugs);
ii = numel(excelsheets_times);

for a = 1:2 %loop metric

    for c = 1:cc %loop cellline
        %     for d = 1:dd %loop drug

        fig = figure;
        fig.Position = [787,627,1401,533];
        sgtitle(celllines{c},'FontSize',14,'FontName','Helvetica Neue','FontWeight','bold');
        %         sgtitle(drug{d},'FontSize',14,'FontName','Helvetica Neue','FontWeight','bold');

        for i = 1:ii %loop input

            %load data
            inputsheet = append(excelsheets_times{i},channel{a});
            [alldata] = readmatrix(file,'sheet',inputsheet);
            %save data from all drugs per cell line
            data(:,i) = alldata(:,c+2);
            %save data from all celllines per cell drug
            %             data(i,:) = alldata(d,:); %

        end

        for k = 1:2 %create two polarhistograms side-by-side

            if k == 1
                subplot(1,2,1) %overlay max and min benefit times with different colors
            else
                subplot(1,2,2) %overlay max and min benefit times + min difference, only show min difference
            end

            for j = 1:3

                if k == 1
                    color = colorall1{j};
                    edgecolor = edgecolor1{j};
                else
                    color = colorall2{j};
                    edgecolor = edgecolor2{j};
                end

                %plot data from all drugs per cell line
                alpha = fliplr(deg2rad(rmmissing(data(:,j)/ 24 * 360)));
                %plot data from all celllines per cell drug
                %                 alpha = deg2rad(rmmissing(data(j,:))/ 24 * 360);

                polarhistogram(alpha,'Normalization','probability','BinEdges',deg2rad(0:30:360),'FaceColor',color,'EdgeColor',edgecolor,'FaceAlpha',0.4); hold all;
                %                 polarhistogram(alpha,'BinEdges',deg2rad(0:30:360),'FaceColor',color,'EdgeColor',edgecolor,'FaceAlpha',0.4); hold all;

                %create legend and fine-tune plot
                if k == 1
                    legend(legendentries,'Location','eastoutside');
                end

                ax = gca;
                ax.PositionConstraint = "innerposition";
                titletext = titles{k};
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
        filename = append('polarhistogram_ToD_times_',celllines{c},'_',channel{a});
        %         filename = append('polarhistogram_ToD_times_',drug{d},'_',metric{a});
        destination = append(pathtofile,'polarhistograms_ToD_times/');
        saveas(fig, [destination filename '.svg']);

    end %loop c cellline

end %loop m metric

end %function