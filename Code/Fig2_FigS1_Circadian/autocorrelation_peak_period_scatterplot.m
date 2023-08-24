function scatterplot_autocorr_peak_period(celllines,file)

%exclude SY5Y because it has periods above 30h

color_string = {'[0,0,0]';'[0.64,0.24,0.98]';'[0.64,0.24,0.98]';'[0.64,0.24,0.98]';'[0.04,0.75,0.1]';'[1.00,0.00,0.00]';'[1.00,0.00,0.00]';'[1.00,0.00,0.00]'};
reporter = {'Bmal1';'Per2'};
font = 'Helvetica Neue';

fig = figure;
fig.Position = [420,285,525,425];

for a = 1:2 %loop reporter

    sheet1 = append(reporter{a},'_lag');
    sheet2 = append(reporter{a},'_peak');

    [lag] = xlsread(file,sheet1);
    [peak] = xlsread(file,sheet2);

    %exclude SY5Y & U2OS Cry1/2-dKO due to periods above the circadain
    %range, exclude U2OS Cry1-sKO from per channel (no data)

    if a == 1
        lag(:,[6,8,9]) = NaN;
        peak(:,[6,8,9]) = NaN;
    elseif a == 2
        lag(:,[6,9]) = NaN;
        peak(:,[6,9]) = NaN;
    end

    x = mean(lag,1,'omitnan');
    x_err = std(lag,[],1,'omitnan');
    y = mean(peak,1,'omitnan');
    y_err = std(peak,[],1,'omitnan');

    ii = numel(color_string);

    for i = 1:ii
        color(i,:) = str2num(color_string{i});
    end

    eb(1) = errorbar(x,y,x_err, 'horizontal', 'LineStyle', 'none'); hold all
    eb(2) = errorbar(x,y,y_err, 'vertical', 'LineStyle', 'none');
    set(eb, 'LineWidth', 0.75, 'Color','black')

    marker = 's^dhv<ox>';
    lgdlocation = 'northeastoutside';

    s = gscatter(x,y,celllines,color,marker); hold all

    if a == 2
        set(s,'MarkerSize',13,{'MarkerFaceColor'},get(s,'Color'));
        legend('Location',lgdlocation,'FontSize',16);
    else
        set(s,'MarkerSize',13,'LineWidth',2,'MarkerFaceColor','w');
    end

    if a == 1
        period_bmal = lag;
        RI_bmal = peak;
    elseif a == 2
        period_per = lag;
        RI_per = peak;
    end

    varstoclear = {'lag';'peak';'x';'x_err';'y';'y_err'};
    clear(varstoclear{:})

end

x = [period_per;period_bmal];
y = [RI_per;RI_bmal];

regr_x = rmmissing(x(:));
regr_y = rmmissing(y(:));

%linear regression
mdl = fitlm(regr_x,regr_y);
plot(mdl,'color','none','HandleVisibility','off','LineWidth',1,'Marker','none');
str =['R^2 = ',sprintf('%.2f',mdl.Rsquared.Ordinary)];
annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black')

hold off

%settings
ylabeltext = 'Autocorrelation 2^{nd} peak';
xlabeltext = 'Lag 2^{nd} peak (h)';

%specify remaining settings and appearance of the plot
ylabel(ylabeltext,'FontSize',16,'FontName',font);
xlabel(xlabeltext,'FontSize',16,'FontName',font);

ax = gca;
box on;
grid on;
set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','LineWidth',0.9,'FontName',font,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

end %function