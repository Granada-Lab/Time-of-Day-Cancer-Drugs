function ranking_circadian_strength_barplot(file,cellline)

%Carolin Ector, 23.08.2023
%%Function ranks cell lines based on their normalized circadian strength values

%Time-of-Day-Cancer-Drugs Manuscript Fig. 2k

%input: stored in ranking_circadian_strength_workspace.mat
% file: excel sheet where detrended circadian time-series data is stored
% celllines: names of the cell lines being analysed

%load data
[num] = xlsread(file,'normalized');

font = 'Helvetica Neue';

meandatanorm = mean(num,1,'omitnan');
stddatanorm = std(num,[],1,'omitnan');

[meansorted, sortIdx] = sort(meandatanorm,'ascend');
datasorted = num(:,sortIdx);
stdsorted = stddatanorm(:,sortIdx);
celllinesorted = cellline(:,sortIdx);

figure;
y = transpose(meansorted);
error = transpose(stdsorted);
bar(y,'FaceColor','flat');
hold on

%adding centered error bars
ngroups = size(y, 1);
nbars = size(y, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), error(:,i), '.','LineWidth', 0.75, 'Color','black');
    s = scatter(x,datasorted,36,'k','filled','MarkerFaceAlpha',0.8','jitter','on','jitterAmount',0.15,'HandleVisibility','off');
    set(s,'MarkerEdgeColor','none');
end

hold off
ax = gca;
set(gca, 'XTickLabel', celllinesorted);
set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','LineWidth',0.9,'FontName',font,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

end %function