function ranking_growth_metrics(file_growthmetrics,cellline)

%Carolin Ector, 25.08.2023

%Function ranks cell lines by their individual growth features (doubling time and growth rate) in bar plots 

%Time-of-Day-Cancer-Drugs Manuscript Fig. 3f

%input: stored in "growth_analysis_workspace.mat"
% file_growthmetrics: excel file where growth data is stored
% cellline: names of the cell lines with growth data

channel = {'CellNr';'Conf'};
channel_long = {'Cell Number';'Confluency (%)'};

for a = 1:2 %loop a channel

    %load data
    sheet = append('Metrics_',channel{a});
    [data] = xlsread(file_growthmetrics,sheet);

    growthrate = data(1:3,:);
    doublingtime = data(6:7,:);

    ylabeltext = 'Doubling time (h)';
    xlabeltext = 'Growth rate (1/h)';

    %% Bar Diagram
    fig = figure;

    x = 1:1:10;

    %sort y-axis by highest to lowest growthrate
    [y1, sortIdx] = sort(growthrate(1,:),'ascend');
    err1_ranked = (growthrate(2,sortIdx));
    err2_ranked = (growthrate(3,sortIdx));
    xaxisnames1 = cellline(sortIdx,:);

    hold on;
    bar(x,y1,'FaceColor', 'r');
    errorbar(x,y1,(y1-err1_ranked),(err2_ranked-y1),'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
    ylabel(xlabeltext);
   
    ax = gca;
    box on;
    grid on;
    set(ax, 'XTick', 1:10, 'XTickLabel', xaxisnames1);
    set(ax,'LineWidth',0.9,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

    filename2 = append('barplot_ranking_growthrate_',channel{a});
%     saveas(fig2, [filename2, '.svg']);

    fig3 = figure;

    %use same sorting as for growthrate to "mirror" the bar plots later on in inkscape
    y2 = (doublingtime(1,sortIdx));
    err3_ranked = (doublingtime(2,sortIdx));
    xaxisnames2 = cellline(sortIdx,:);

    hold on;
    bar(x,y2,'FaceColor', 'b');
    errorbar(x,y2,err3_ranked,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
    ylabel(ylabeltext);
   
    ax = gca;
    box on;
    grid on;
    set(ax, 'XTick', 1:10, 'XTickLabel', xaxisnames2);
    set(ax,'LineWidth',0.9,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

    filename3 = append('barplot_ranking_doublingtime_',channel{a});
%     saveas(fig3, [filename3, '.svg']);

end %loop a channel
end %function


