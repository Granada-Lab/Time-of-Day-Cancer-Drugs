function correlation_ToDMR_smoothed_unsmoothed(drugs,ToDMR_raw,ToDMR_spline)

%Carolin Ector, 23.08.2023
%calculate the linear correlation coefficient between ToDMR values (maximum ToD response range) of unsmoothed and smoothed reponse data 
%Time-of-Day-Cancer-Drugs Manuscript Fig. S7

%input: stored in ToD_analysis_workspace.mat
% drugs: names of drugs the cell lines have been treated with
% ToD_raw: ToDMR values from measured datapoints (unsmoothed)
% ToD_spline: ToDMR values from spline-smoothed datapoints

font = 'Helvetica Neue';
color = hsv(9);

%exclude HCC1937_2 measurements
ToDMR_raw(:,[1,9]) = [];
ToDMR_spline(:,[1,9]) = [];

x = cell2mat(table2array(ToDMR_spline));
y = table2array(ToDMR_raw);

fig = figure;
fig.Position = [420,285,525,425];

xlabeltext = 'ToD Max. Range Spline Smooth';
ylabeltext = 'ToD Max. Range Datapoints';
 
for d = 1:10
    s = gscatter(x(:,d),y(:,d),drugs,color); hold all
    set(s,'MarkerSize',19,{'MarkerFaceColor'},get(s,'Color'));
end

legend('Location','northeastoutside','FontSize',16);

%linear regression
regr_x = rmmissing(x(:));
regr_y = rmmissing(y(:));
mdl = fitlm(regr_x,regr_y);
plot(mdl,'color','none','HandleVisibility','off','LineWidth',1,'Marker','none');
str =['R^2 = ',sprintf('%.2f',mdl.Rsquared.Ordinary)];
annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black')

hold off

ylim([-0.04 1.24]);
xlim([-0.04 1.24]);
xticks (0:0.3:1.2);

xlabel(xlabeltext,'FontSize',16,'FontName',font);
ylabel(ylabeltext,'FontSize',16,'FontName',font);

ax = gca;
box on;
grid on;
set(ax,'LineWidth',0.9,'FontName',font,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

%save figure
saveas(fig, [filename, '.svg']);

end %function