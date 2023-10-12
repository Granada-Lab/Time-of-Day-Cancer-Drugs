function total_dominance_ranking_boxplot

%Carolin Ector, 23.08.2023

% Function ranks overall contribution (total dominance) of each metric to the ToDMR of all drugs combined

% input data table is stored in corresponding Figshare folder

%Time-of-Day-Cancer-Drugs Manuscript Fig. 5i

inputfile = 'Data_Fig5/TotalDominanceMRSpline_Bmal1_plusgrowthrate_Data.xlsx';

fig = figure;
fig.Position = [420,285,525,425];

[data,text,~] = xlsread(inputfile);

values = text(1,2:end);
data = data*100; %percentage of contribution

means = mean(data,1);

% med = median(data);
[meansAscend, sortIdx] = sort(means,'descend');
dataAscend = data(:,sortIdx);
valuesAscend = values(:,sortIdx);

boxplot(dataAscend,'Labels',valuesAscend); hold all

datasize = size(dataAscend);
x = repmat(1:datasize(:,2),datasize(:,1),1);
s = scatter(x,dataAscend,60,'k','filled','MarkerFaceAlpha',0.8','jitter','on','jitterAmount',0.15);
set(s,'MarkerEdgeColor',[0 0 0]);

% Add mean markers to the box plot
scatter(x, meansAscend, 60, 'r', 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
hold off;

ylabel('Contribution (%)')

ax = gca;
set(findobj(gca,'type','line'),'linew',1);
box on;
grid on
set(ax,'YLimitMethod', 'padded','LineWidth',0.9,'FontName','Helvetica Neue','FontSize',18,'XMinorGrid','off','YMinorGrid','off');

end %function
