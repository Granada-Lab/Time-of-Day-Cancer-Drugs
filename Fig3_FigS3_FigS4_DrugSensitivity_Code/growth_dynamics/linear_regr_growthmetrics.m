function linear_regr_growthmetrics(file_growthmetrics,regression_combi,cellline,color_all)

%Carolin Ector, 25.08.2023

%Function correlates diferent growth metrics by linear regression analysis
%1. doubling times obtained from the fluorescent channel (cell number)
% and the brieghtfield channel (confluency)
%2. doubling times and growth rates from exponential curve fit

%Time-of-Day-Cancer-Drugs Manuscript Fig.S2b, Fid.S2c

%input: stored in growth_analysis_workspace.mat
% file: excel file where growth data is stored
% cellline: names of the cell lines with growth data
% color_all: color array containing a distinct color for each cell line

cc = numel(cellline);

%load data
[cellnr] = xlsread(file_growthmetrics,'Metrics_CellNr');
[conf] = xlsread(file_growthmetrics,'Metrics_Conf');

%convert color_string to num (scatterplot can't read color_string)
for i = 1:cc
    color(i,:) = str2num(color_all{i,:});
end

for b = 1:2

    combination = regression_combi{b};

    fig = figure;
    fig.Position = [420,285,525,425];

    %linear regression analysis
    hold all
    if b == 1 %doublingtime cellnr vs conf
        regr_x = rmmissing(cellnr(6,:));
        regr_y = rmmissing(conf(6,:));
        combinationtext = 'doublingtime_cellnr_vs_conf';
    elseif b == 2 %growthrate cellnr vs doublingtime cellnr
        regr_x = rmmissing(cellnr(1,:));
        regr_y = rmmissing(cellnr(6,:));
        combinationtext = 'growthrate_vs_doublingtime_cellnr';
    end
    mdl = fitlm(regr_x,regr_y);
    plot(mdl,'color','none','HandleVisibility','off','LineWidth',1,'Marker','none');

    coeffs = mdl.Coefficients;
    pvalue = table2array(coeffs(2,4));

    str =['R^2 = ',sprintf('%.2f',mdl.Rsquared.Ordinary),'  p = ',sprintf('%.4f',pvalue)];
    annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','none','FontSize',16);

    if b == 1 %doublingtime cellnr vs conf
        eb(1) = errorbar(cellnr(6,:),conf(6,:),cellnr(7,:),'horizontal', 'LineStyle', 'none', 'HandleVisibility','off');
        eb(2) = errorbar(cellnr(6,:),conf(6,:),conf(7,:),'vertical','LineStyle','none','HandleVisibility','off');
        s = gscatter(cellnr(6,:),conf(6,:),cellline,color);
    elseif b == 2 %growthrate cellnr vs doublingtime cellnr
        eb(1) = errorbar(cellnr(1,:),cellnr(6,:),(cellnr(1,:)-cellnr(2,:)),(cellnr(3,:)-cellnr(1,:)),'horizontal', 'LineStyle', 'none', 'HandleVisibility','off');
        eb(2) = errorbar(cellnr(1,:),cellnr(6,:),cellnr(7,:),'vertical','LineStyle','none','HandleVisibility','off');
        s = gscatter(cellnr(1,:),cellnr(6,:),cellline,color);
    end

    set(eb, 'LineWidth', 0.75, 'Color','k');

    set(s,'MarkerSize',30,{'MarkerFaceColor'},get(s,'Color'));
    legend('off'); % Disable legend for the mdl plot

    hold off

    %specify remaining settings and appearance of the plot
    xlabel(combination{1});
    ylabel(combination{2});

    ax = gca;
    set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded');

    box on;
    grid on;
    set(ax,'LineWidth',0.9,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

    %save figure
    filename = append('linear_regression_',combinationtext);
    % saveas(fig, [filename, '.svg']);

end %function