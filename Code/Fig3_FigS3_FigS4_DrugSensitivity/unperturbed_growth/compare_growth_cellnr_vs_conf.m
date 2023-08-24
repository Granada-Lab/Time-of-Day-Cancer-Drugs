function compare_growth_cellnr_vs_conf(celllineselected,colorselected)

cc = numel(celllineselected);

%convert color_string to num (scatterplot can't read color_string)
for i = 1:cc
    color(i,:) = str2num(colorselected{i,:});
end

    %load data
    folder = append('selected_celllines_manuscript/');
    inputfile1 = append(folder,'Results_Growth_Data_48_well_Selection_Manuscript_CellNr_short.xlsx');
    inputfile2 = append(folder,'Results_Growth_Data_48_well_Selection_Manuscript_Conf_short.xlsx');

    [cellnr] = xlsread(inputfile1,'doublingtime');
    [conf] = xlsread(inputfile2,'doublingtime');


    ylabeltext = 'Confluence';
    xlabeltext = 'Cell Number';
    titletext = 'Doubling Time [hours]';
%     titletext = 'Growth Rate [1/h]';

    fig = figure;
    fig.Position = [420,285,525,425];

    %linear regression
    hold all
    regr_x = rmmissing(cellnr(1,:));
    regr_y = rmmissing(conf(1,:));
    mdl = fitlm(regr_x,regr_y);
    plot(mdl,'color','none','HandleVisibility','off','LineWidth',1,'Marker','none');
    
    coeffs = mdl.Coefficients;
    pvalue = table2array(coeffs(2,4));

    str =['R^2 = ',sprintf('%.2f',mdl.Rsquared.Ordinary),'  p = ',sprintf('%.4f',pvalue)];
    annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','none','FontSize',16);

    %doublingtime
    eb(1) = errorbar(cellnr(1,:),conf(1,:),cellnr(2,:),'horizontal', 'LineStyle', 'none', 'HandleVisibility','off');
    eb(2) = errorbar(cellnr(1,:),conf(1,:),conf(2,:),'vertical','LineStyle','none','HandleVisibility','off');
    
    %growthrate
%     eb(1) = errorbar(cellnr(1,:),conf(1,:),(cellnr(1,:)-cellnr(2,:)),(cellnr(3,:)-cellnr(1,:)),'horizontal', 'LineStyle', 'none', 'HandleVisibility','off');
%     eb(2) = errorbar(cellnr(1,:),conf(1,:),(conf(1,:)-conf(2,:)),(conf(3,:)-conf(1,:)),'vertical','LineStyle','none','HandleVisibility','off');
    set(eb, 'LineWidth', 0.75, 'Color','k');

    s = gscatter(cellnr(1,:),conf(1,:),celllineselected,color);
    set(s,'MarkerSize',30,{'MarkerFaceColor'},get(s,'Color'));
    legend('off'); % Disable legend for the mdl plot
%     legend(celllineselected,'Location','northeastoutside','FontSize',16);

    hold off

    %specify remaining settings and appearance of the plot
    title(titletext,'FontWeight','normal');
    ylabel(ylabeltext);
    xlabel(xlabeltext);

    ax = gca;
    set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded');

    box on;
    grid on;
    set(ax,'LineWidth',0.9,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

    %save figure
    filename = append(folder,'scatterplot_doublingtime_cellnr_vs_conf');
    saveas(fig, [filename, '.svg']);

end %function


