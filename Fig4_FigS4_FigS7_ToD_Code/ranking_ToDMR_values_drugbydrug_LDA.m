function ranking_ToDMR_values_drugbydrug_LDA(file,channel,celllines_LDA,drugs)

%Carolin Ector, 23.08.2023

%Time-of-Day-Cancer-Drugs Manuscript Fig. 6e & S6b-h

%rank cell lines by their ToDMR values per drug to visualize the grouping into high and low ToD sensitivity groups by LDA (bar colors adjusted later in inkscape)

%input: stored in "ToD_analysis_workspace.mat"
% file: excel file where ToDMR values are stored
% celllines_LDA: names of the cell lines being analysed for LDA (differently ordered than excel sheets for max and min benefit foldchanges/ToD)
% channel: either brightfield (confluency) or fluorescent (cell number) channel from live-imaging, in the manuscript we focus on the fluorescent channel 
% drugs: names of drugs the cell lines have been treated with

dd = numel(drugs);

%% Load and normalize data
for a = 1:2 %loop a channel

    sheet = append('ToDMR_',channel{a});

    [alldata] = readmatrix(file,'Sheet',sheet);
    alldata([9,6],:) = []; %exclude MCF10A and HCC1937_2
    alldata(:,[1,end]) = []; %exclude first row (all NaN) and Olaparib (only 2 celllines)

    for d = 1:dd-1 %loop d drugs (exclude Olaparib)

        data = alldata(:,d);
        ylabeltext = celllines_LDA;

        if d == 3 %no data for SUM149PT + Alpelisib
            data(9,:) = [];
            ylabeltext(9,:) = [];
        elseif d==8 %no data for HCC1806 + Cisplatin
            data(4,:) = [];
            ylabeltext(4,:) = [];
        end

        x = (1:1:numel(ylabeltext))';
        [ybar, sortIdx] = sort(data,'descend');
        celllinesorted = ylabeltext(sortIdx,:);
        
        median_data = median(data);

        %add median (=cutoff for binarization) to the title
        titletext = append(drugs{d},' M = ',num2str(median_data));

        fig = figure;
        fig.Position = [1189,918,261,420];
        barh(x,ybar,'BarWidth',0.7,'FaceColor',[0.9,0.9,0.9],'LineWidth',0.9); hold on

        yticklabels(celllinesorted);
        title(titletext);
        ax = gca;
        ylim([0.4,(numel(ylabeltext)+0.6)]);
        axOpt = {'linewidth',1.5,'YGrid','off','XGrid','off','Box','on','Color','none','FontSize',16};
        set(gca,'XLimitMethod', 'padded');
        set(ax,axOpt{:});

        filename = append('barplot_ToDMR_',drugs{d},'_sortedbymedian_LDA');
        saveas(fig, filename, 'svg');

        clear data
        clear ylabeltext
    end

end %loop output
end %function
