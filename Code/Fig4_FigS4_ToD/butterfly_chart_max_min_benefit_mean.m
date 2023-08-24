function butterfly_chart_max_min_benefit_mean(file,excelsheets,celllines,channel,drugs)

%Carolin Ector, 23.08.2023
%Time-of-Day-Cancer-Drugs Manuscript Fig. 4k

%create butterfly charts to evaluate the extend of the difference in responses between the cancer cell line and MCF10A 
%at the respective maximum and minimum ToD treatment times --> merged for all cell lines per drug or vice versa
% - maximum benefit (highest response difference between cancer cell line and mcf10a with higher toxicity against cancer cell line as opposed to mcf10a)
% - minimum benefit (highest response difference between cancer cell line and mcf10a with higher toxicity against mcf10a as opposed to cancer cell line)

%input: stored in ToD_analysis_workspace.mat
% file: excel file where maximum and minumum benefit foldchanges and times are stored
% excelsheets: sheet names for foldchanges in responses at maximum and minimum benefit
% celllines: names of the cell lines being analysed
% channel: either brightfield (confluency) or fluorescent (cell number) channel from live-imaging, in the manuscript we focus on the fluorescent channel 
% drugs: names of drugs the cell lines have been treated with

ii = numel(excelsheets);
celllines(end-3,:) = []; %remove HCC1937_2
drugs(end,:) = []; %remove Olaparib (only 1 cancer cell line tested)

for a = 1:2 %loop a channel

    for i = 1:ii %loop i input

        %load data
        inputsheet = append(excelsheets{i},channel{a});
        [alldata] = readmatrix(file,'sheet',inputsheet);

        %remove HCC1937_2
        alldata(:,end-3) = []; 

        %average data
        % by drug
        meandrug(:,i) = mean((1-alldata),2,'omitnan');
        stddrug(:,i) = std((1-alldata),[],2,'omitnan');
        % by cellline
        meancell(:,i) = rmmissing(mean((1-alldata),1,'omitnan'));
        stdcell(:,i) = rmmissing(std((1-alldata),[],1,'omitnan'));

    end %loop i input

    %remove Olaparib (only 1 cancer cell line tested)
    meandrug(9,:) = [];
    stddrug(9,:) = [];

    for k = 1:2 %loop k different plots (mean drug or celllines)

        if k == 1 %mean drug
            x = transpose((1:1:numel(drugs)));
            y = meandrug;
            yaxisnames = drugs;
            name = 'drugs';
            err = stddrug;
        else %mean celllines
            x = transpose((1:1:numel(celllines)));
            y = meancell;
            yaxisnames = celllines;
            name = 'celllines';
            err = stdcell;
        end

        fig = figure;

        %sort y-axis by highest to lowest maximum benefit
        [y2, sortIdx] = sort(y(:,2),'ascend');
        y1 = y(:,1);
        y1 = y1(sortIdx,:);
        err = err(sortIdx,:);
        yaxisnames = yaxisnames(sortIdx,:);

        % Create the butterfly chart
        hold all; 
        barh(x, y2, 'FaceColor', [0.549 0.9569 0.7373], 'BarWidth', 0.6,'LineWidth',0.8); %First set of bars
        barh(x, y1, 'FaceColor', [0.9843 0.5804 0.5804], 'BarWidth', 0.6,'LineWidth',0.8);  % Second set of bars (negative values)
        
        % Finding the number of groups and the number of bars in each group
        ngroups = size(y1, 1);
        nbars = size(y1, 2);
        % Calculating the width for each bar group
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        % Set the position of each error bar in the centre of the main bar
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        for i = 1:nbars
            % Calculate center of each bar
            x2 = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(y2,x2,err(:,2),err(:,2),'horizontal','k', 'linestyle', 'none','LineWidth',0.8);
            errorbar(y1,x2,err(:,1),err(:,1),'horizontal','k', 'linestyle', 'none','LineWidth',0.8);
        end

        hold off;

        % Adjust axis limits and labels
        ylimits = size(y,1)+0.5;
        ylim([0.5, ylimits]);
        yticks(1:1:ylimits);
        yticklabels(yaxisnames);

        xlabeltext = append('Mean Foldchange ToDMR');
        xlabel(xlabeltext);

        % Title and legend
        title('Chronotherapeutic Effect Size','FontSize',22);
        lgd = legend('Max Benefit', 'Min Benefit');
        set(lgd,'Location','southeast');

        ax = gca;
        set(ax,'XLimitMethod','padded','linewidth',1.5,'Color','w','FontSize',22,'FontName','Helvetica Neue','box','on');

        %save figure
        filename = append('butterfly_max_min_ToD_effect_',name,'_mean_',channel{a});
        saveas(fig, [destination filename '.svg']);

        valuestoclear = {'data';'drug';'sortIdx';'x';'y1';'y2';'err'};
        clear(valuestoclear{:});

    end %loop k different plots

end %loop a channel
end %function