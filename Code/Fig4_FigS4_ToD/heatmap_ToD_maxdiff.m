function heatmap_ToD_maxdiff(file,celllines_heatmap,channel,drugs,color)

%Carolin Ector, 23.08.2023

%Time-of-Day-Cancer-Drugs Manuscript Fig. 4f+g

%create clustergram of ToDMR values and create adjacent barplots showing a rank of ToDMR values per drug or per cellline 

%input: stored in ToD_analysis_workspace.mat
% file: excel file where ToDMR values are stored
% celllines_heatmap: names of the cell lines included in the heatmap and ranking
% channel: either brightfield (confluency) or fluorescent (cell number) channel from live-imaging, in the manuscript we focus on the fluorescent channel 
% drugs: names of drugs the cell lines have been treated with
% color: color for clustergram map

drugs(end,:) = [];
ybarmax_d = 8.6;
ybarmax_c = 10.6;

% Load and normalize data
for a = 1:2 %loop a channel

    sheet = append('ToDMR_',channel{a});

    [data] = readmatrix(file,'Sheet',sheet);
    data(9,:) = []; %exclude HCC1937_2
    data(:,[1,end]) = []; %exclude first row (all NaN) and Olaparib (only 2 celllines)

    %calculate mean ToDMR values by drug or by cellline
    mean_cell = mean(data,2,'omitnan');
    std_cell = std(data,[],2,'omitnan');
    mean_drug = mean(data,1,'omitnan');
    std_drug = std(data,[],1,'omitnan');

    %sort data from max to min ToDMR for adjacent barplots
    [ybar1, sortIdx_d] = sort(mean_drug,'descend');
    [ybar2, sortIdx_c] = sort(mean_cell,'descend');

    ebar1 = std_drug(:,sortIdx_d);
    ebar2 = std_cell(sortIdx_c,:);
    drugsorted = drugs(sortIdx_d,:);
    celllinesorted = celllines_heatmap(sortIdx_c,:);

    %% clustergram of ToDMR values

    Y=data';
    cgo = clustergram(Y,'Colormap',color,'ImputeFun', @knnimpute);
    set(cgo,'RowLabels',drugs,'ColumnLabels',celllines_heatmap,'AnnotColor','k','AnnotPrecision',1);
    set(cgo,'ColumnLabelsRotate',90,'RowLabelsRotate',0,'Displayrange',0.45);

    %% ranking of drugs by their average ToDMR value
    fig = figure;

    ybar_d = flip(ybar1);
    ebar_d = flip(ebar1);
    xlabels = flip(drugsorted);
    n = numel(ebar_d);
    ebarlow(1:n,:) = 0;

    x_d = 1:1:numel(drugsorted);

    bar(x_d,ybar_d,'BarWidth',0.7,'FaceColor',[0.9,0.9,0.9],'LineWidth',0.9); hold on
    xticklabels(xlabels);
    ax = gca;
    ylim([0.4,ybarmax_d]);
    axOpt = {'linewidth',0.8,'YGrid','off','XGrid','off','Box','on','Color','none','FontSize',14};
    set(gca,'YTickLabel',[],'XLimitMethod', 'padded','InnerPosition',[0.41,0.11,0.07,0.81],'OuterPosition',[0.40,0,0.11,0.99]);
    set(ax,axOpt{:});

    % Finding the number of groups and the number of bars in each group
    ngroups = size(ybar_d, 2);
    nbars = size(ybar_d, 1);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*1-1) * groupwidth / (2*nbars);
    errorbar(ybar_d,x,ebarlow,ebar_d,'horizontal','k', 'linestyle', 'none','LineWidth',0.8);

    hold off

    %% ranking of celllines by their average ToDMR value
    figure;

    ybar_c = flip(ybar2);
    ebar_c = flip(ebar2);
    n2 = numel(ebar_c);
    ebarlow2(1:n2,:) = 0;

    x_c = (1:1:numel(celllinesorted))';

    bar(x_c,ybar_c,'BarWidth',0.7,'FaceColor',[0.9,0.9,0.9],'LineWidth',0.9); hold on
    xticklabels(celllinesorted);
    ax = gca;
    ylim([0.4,ybarmax_c]);
    axOpt = {'linewidth',0.8,'YGrid','off','XGrid','off','Box','on','Color','none','FontSize',14};
    set(gca,'YTickLabel',[],'XLimitMethod', 'padded','InnerPosition',[0.410797101449275,0.113809523220227,0.077714387062213,0.811190476779773], ...
        'OuterPosition',[0.589729367425476,0,0.113317587080332,0.999649430378508]);
    set(ax,axOpt{:});

    % Finding the number of groups and the number of bars in each group
    ngroups = size(ybar_c, 2);
    nbars = size(ybar_c, 1);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*1-1) * groupwidth / (2*nbars);
    errorbar(ybar_c,x,ebarlow2,ebar_c,'horizontal','k', 'linestyle', 'none','LineWidth',0.8);

    hold off

    %     filename = append('heatmap_ToD_cosamplitude_allcelllines_alldrugs_',metric{a});
    %     saveas(fig, filename, 'svg');
    %     saveas(fig, filename, 'fig');
    %
    %     clear alldata

end %loop a channel

end %function
