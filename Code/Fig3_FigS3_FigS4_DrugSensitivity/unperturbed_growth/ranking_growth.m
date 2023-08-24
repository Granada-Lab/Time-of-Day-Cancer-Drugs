function ranking_growth(celllineselected,colorselected)

metric = {'CellNr';'Conf'};
metriclong = {'Cell Number';'Confluence'};
cc = numel(celllineselected);

%convert color_string to num (scatterplot can't read color_string)
for i = 1:cc
    color(i,:) = str2num(colorselected{i,:});
end

for m = 1:2

    %load data
    folder = append('selected_celllines_manuscript/');
    inputfile = append(folder,'Results_Growth_Data_48_well_Selection_Manuscript_',metric{m},'_short.xlsx');

    [growthrate] = xlsread(inputfile,'growthrate');
    [doublingtime] = xlsread(inputfile,'doublingtime');

    ylabeltext = 'Doubling Time [hours]';
    xlabeltext = 'Growth Rate [1/h]';
    titletext = metriclong{m};

    %% Scatterplot
%     fig = figure;
%     fig.Position = [420,285,525,425];
% 
%     %linear regression
%     hold all
%     regr_x = rmmissing(growthrate(1,:));
%     regr_y = rmmissing(doublingtime(1,:));
%     mdl = fitlm(regr_x,regr_y);
%     plot(mdl,'color','none','HandleVisibility','off','LineWidth',1,'Marker','none');
% 
%     coeffs = mdl.Coefficients;
%     pvalue = table2array(coeffs(2,4));
% 
%     str =['R^2 = ',sprintf('%.2f',mdl.Rsquared.Ordinary),'  p = ',sprintf('%.4f',pvalue)];
%     annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black')
% 
%     eb(1) = errorbar(growthrate(1,:),doublingtime(1,:),(growthrate(1,:)-growthrate(2,:)),(growthrate(3,:)-growthrate(1,:)),'horizontal', 'LineStyle', 'none', 'HandleVisibility','off');
%     eb(2) = errorbar(growthrate(1,:),doublingtime(1,:),doublingtime(2,:),'vertical','LineStyle','none','HandleVisibility','off');
%     set(eb, 'LineWidth', 0.75, 'Color','k');
% 
%     s = gscatter(growthrate(1,:),doublingtime(1,:),celllineselected,color);
%     set(s,'MarkerSize',30,{'MarkerFaceColor'},get(s,'Color'));
%     legend('off'); % Disable legend for the mdl plot
%     %     legend(celllineselected,'Location','northeastoutside','FontSize',16);
% 
%     hold off
% 
%     %specify remaining settings and appearance of the plot
%     title(titletext,'FontWeight','normal');
%     ylabel(ylabeltext);
%     xlabel(xlabeltext);
% 
%     ax = gca;
%     set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded');
% 
%     box on;
%     grid on;
%     set(ax,'LineWidth',0.9,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');
% 
%     %save figure
%     filename = append(folder,'scatterplot_growthrate_vs_doublingtime_',metric{m});
%     %     saveas(fig, [filename, '.svg']);

    %% Bar Diagram
    fig2 = figure;

    x = 1:1:10;

    %sort y-axis by highest to lowest growthrate
    [y1, sortIdx] = sort(growthrate(1,:),'ascend');
    err1_ranked = (growthrate(2,sortIdx));
    err2_ranked = (growthrate(3,sortIdx));
    xaxisnames1 = celllineselected(sortIdx,:);

    hold on;
    bar(x,y1,'FaceColor', 'r');
    errorbar(x,y1,(y1-err1_ranked),(err2_ranked-y1),'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
    ylabel(xlabeltext);
   
    ax = gca;
    box on;
    grid on;
    set(ax, 'XTick', 1:10, 'XTickLabel', xaxisnames1);
    set(ax,'LineWidth',0.9,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

    filename2 = append(folder,'barplot_ranking_growthrate_',metric{m});
%     saveas(fig2, [filename2, '.svg']);

    fig3 = figure;

%     %sort y-axis by highest to lowest doublingtime
%     [y2, sortIdx] = sort(doublingtime(1,:),'ascend');
%     err3_ranked = (doublingtime(2,sortIdx));
%     xaxisnames2 = celllineselected(sortIdx,:);

    y2 = (doublingtime(1,sortIdx));
    err3_ranked = (doublingtime(2,sortIdx));
    xaxisnames2 = celllineselected(sortIdx,:);

    hold on;
    bar(x,y2,'FaceColor', 'b');
    errorbar(x,y2,err3_ranked,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
    ylabel(ylabeltext);
   
    ax = gca;
    box on;
    grid on;
    set(ax, 'XTick', 1:10, 'XTickLabel', xaxisnames2);
    set(ax,'LineWidth',0.9,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

    filename3 = append(folder,'barplot_ranking_doublingtime_',metric{m});
%     saveas(fig3, [filename3, '.svg']);

    %% Butterfly Chart

%     fig4 = figure;
% 
%     x = 1:1:10;
% 
%      %sort y-axis by highest to lowest growthrate
%     [y1, sortIdx] = sort(growthrate(1,:),'ascend');
%     err1_ranked = growthrate(2,sortIdx);
%     err2_ranked = growthrate(3,sortIdx);
%     y2 = doublingtime(1,sortIdx);
%     err3_ranked = doublingtime(2,sortIdx);
%     xaxisnames1 = celllineselected(sortIdx,:);
% 
%     hold all;
%     barh(x, y2, 'FaceColor', [0.549 0.9569 0.7373], 'BarWidth', 0.6,'LineWidth',0.8); %First set of bars
%     barh(x, -y1, 'FaceColor', [0.9843 0.5804 0.5804], 'BarWidth', 0.6,'LineWidth',0.8);  % Second set of bars (negative values)
% 
%     % Finding the number of groups and the number of bars in each group
%     ngroups = size(y1, 1);
%     nbars = size(y1, 2);
%     % Calculating the width for each bar group
%     groupwidth = min(0.8, nbars/(nbars + 1.5));
%     % Set the position of each error bar in the centre of the main bar
%     % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
%     for i = 1:nbars
%         % Calculate center of each bar
%         x2 = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%         errorbar(y2,x,err3_ranked,'horizontal','k', 'linestyle', 'none','LineWidth',0.8);
%         errorbar(-y1,x,-(y1-err1_ranked),-(err2_ranked-y1),'horizontal','k', 'linestyle', 'none','LineWidth',0.8);
%     end
% 
%     hold off;


end % metric

end %function


