function total_dominance_ranking_boxplot

%Carolin Ector, 23.08.2023
% rank overall contribution (total dominance) of each metric to the ToDMR of all drug combined

path = 'Total_Dominance_Selected_Features.xlsx';

fig = figure;
fig.Position = [420,285,525,425];

[data,text,~] = xlsread(path);

drug = text(2:end,1);
values = text(1,2:end);
data = data*100; %percentage of contribution

med = median(data,'omitnan');
[~, sortIdx] = sort(med,'ascend');
dataAscend = data(:,sortIdx);
valuesAscend = values(:,sortIdx);

boxplot(dataAscend,'Labels',valuesAscend); hold on

datasize = size(dataAscend);
x = repmat(1:datasize(:,2),datasize(:,1),1);
s = scatter(x,dataAscend,60,'k','filled','MarkerFaceAlpha',0.8','jitter','on','jitterAmount',0.15);
set(s,'MarkerEdgeColor',[0 0 0]);

ylabel('Contribution (%)')

ax = gca;
set(findobj(gca,'type','line'),'linew',1);
box on;
grid on
set(ax,'YLimitMethod', 'padded','LineWidth',0.9,'FontName','Helvetica Neue','FontSize',18,'XMinorGrid','off','YMinorGrid','off');

end %function
