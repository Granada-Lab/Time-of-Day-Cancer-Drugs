function overall_correlation_ranking_barplot(labels,mean,stderr)

%Carolin Ector, 23.08.2023

% rank overall correlations between each metric and the ToD maximum range in sensitivity (ToDMR)
% rank either by metric or by drug

%input: stored in overall_linear_correlation_ranking_workspace.mat
% labels = labels for x-axis
% mean = mean absolute linear correlation values per drug or metric
% stderr = standard error of absolute linear correlation values per drug or metric

for i = 1:2 %loop i values for ranking (either by drug or metric)

    %load data, i=1: by drug, i=2: by metric
    labelnames = labels{i};
    meandata = mean{i};
    stddata = stderr{i};

    fig = figure;

    x = 1:1:(numel(meandata));

    %sort y-axis by highest to lowest absolute correlation
    [y, sortIdx] = sort(meandata,'ascend');
    err_ranked = (stddata(sortIdx,:));
    err_0(1:numel(meandata),1) = 0;
    label = labelnames(sortIdx,:);

    hold on;
    bar(x,y,'FaceColor', 'r');
    errorbar(x,y,err_0,err_ranked,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
    ylabel('Absolute Correlation');

    ax = gca;
    box on;
    grid on;
    ylim([0.2 Inf]);
    set(ax, 'XTick', 1:1:numel(meandata), 'XTickLabel', label);
    set(ax,'LineWidth',1.5,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

end %loop i values for ranking
end %function
