function [coeffs] = correlation_ToDMR_smoothed_unsmoothed(drug,color_drugs)

%Carolin Ector, 23.08.2023

%Function calculates the linear correlation coefficient between ToDMR values (maximum ToD response range) of unsmoothed and smoothed reponse data 

%Time-of-Day-Cancer-Drugs Manuscript Fig. S7

%input: stored in "ToD_overlay_workspace.mat"
% drugs: names of drugs the cell lines have been treated with
% color_drugs: color-coding for drugs used in the manuscript

font = 'Helvetica Neue';
color = color_drugs;
channel = 'CellNr';
dd = numel(drug);
path = '/Users/carolinector/Nextcloud/Manuscripts/ToDMethods/code_and_data_to_submit/Data/Fig4_FigS4_ToD_Data/ToD_experiments_results/';

%load ToDMR values calculated from raw ToD responses  
for d = 1:dd
    sheet_raw = append(drug{d},'_',channel);
    inputfile_raw = append(path,'ToD_maxrange_sorted.xlsx');
    [raw] = readmatrix(inputfile_raw,'Sheet',sheet_raw);
    if d == dd
        raw_sorted(dd,2) = raw(3,2);
        raw_sorted(dd,6) = raw(7,2);
        raw_sorted(dd,[1,3:5,7:11]) = NaN;
    else
        raw_sorted(d,:) = raw(:,2)';
    end
end

raw_sorted(:,9) = []; %exclude HCC1937_2

%load ToDMR values calculated from spline smoothed ToD responses  
sheet_spline = append('ToDMR_spline_',channel);
inputfile_spline = append(path,'ToD_values_splinesmooth_benefit.xlsx');
[spline] = readmatrix(inputfile_spline,'Sheet',sheet_spline);
spline(:,1) = []; %exclude first row (all NaNs)
spline(9,:) = []; %exclude HCC1937_2

x = spline';
y = raw_sorted;

fig = figure;
fig.Position = [420,285,525,425];

xlabeltext = 'ToDMR Spline Smooth';
ylabeltext = 'ToDMR Datapoints';
 
for d = 1:10
    s = gscatter(x(:,d),y(:,d),drug,color); hold all
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

ylim([-0.03 0.93]);
xlim([-0.03 0.93]);
xticks (0:0.3:1.2);
yticks (0:0.3:1.2);

xlabel(xlabeltext,'FontSize',16,'FontName',font);
ylabel(ylabeltext,'FontSize',16,'FontName',font);

ax = gca;
box on;
grid on;
set(ax,'LineWidth',0.9,'FontName',font,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

coeffs = mdl.Coefficients;

%save figure
% saveas(fig, [filename, '.svg']);

end %function