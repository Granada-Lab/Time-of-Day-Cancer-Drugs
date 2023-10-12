function mra_scatterplot(celllines,file,row_bmal,row_per,rep_bmal,rep_per)

%Carolin Ector, 23.08.2023

%Function plots multi-resolution analysis (MRA) values in a scatterplots
%(noise vs. circadian component)

%Time-of-Day-Cancer-Drugs Manuscript Fig. 2j

%input: stored in "workspace_mra_scatterplot.mat"
% file: excel file where MRA data is stored
% celllines: names of the cell lines being analysed
% row_bmal / row_per: row of first replicate of respective cellline for each reporter
% rep_bmal / rep_per: number of replicates per reporter cell line

%assign each cell line a color based on their tissue type
color_string = {'[0,0,0]';'[0.64,0.24,0.98]';'[0.64,0.24,0.98]';'[0.64,0.24,0.98]';'[0.04,0.75,0.11]';'[0.04,0.75,0.11]';'[1.00,0.00,0.00]';'[1.00,0.00,0.00]';'[1.00,0.00,0.00]'};
font = 'Helvetica Neue';

%convert color_string to num (scatterplot can't read color_string)
for i = 1:numel(color_string)
    color(i,:) = str2num(color_string{i,:});
end

%load MRA values obtained from python script: mra_extract_components.py
[num] = xlsread(file);

for a = 1:2 %loop a reporter

    if a == 1
        cc = numel(celllines);
    elseif a == 2
        cc = numel(celllines(1:7,:));
    end

    for c = 1:cc %loop cellline

        if a == 1
            rr = rep_bmal(:,c);
        else
            rr = rep_per(:,c);
        end

        for r = 1:rr %loop r replicate
            if a == 1 %bmal
                row = row_bmal(:,c)+(r-1);
                saverow = 1;
            else %per
                row = row_per(:,c)+(r-1);
                saverow = 2;
            end
            noisevalues(:,r) = num(row,1);
            circvalues(:,r) = num(row,3);
        end

        %create array of mean and std values: columns = celllines, row1 = Bmal1_mean, row2=Per2_mean,
        noise(saverow,c) = mean(noisevalues,2);
        noise_std(saverow,c) = std(noisevalues,[],2);
        circ(saverow,c) = mean(circvalues,2);
        circ_std(saverow,c) = std(circvalues,[],2);

        if rr < 9
            circvalues(:,end+1:9) = NaN;
        end

        if a == 1 && c < 8
            allcirc_bmal(:,c) = transpose(circvalues);
        elseif a == 2 && c < 8
            allcirc_per(:,c) = transpose(circvalues);
        end

        %calculate ratios
        %columns = celllines, row1 = ratio circ/noise bmal, row2 = ratio circ/noise per
        if c < 8
            ratio_cn(saverow,c) = circ(saverow,c)/noise(saverow,c);
            %ratio_cn_std(saverow,c) = noise_std(saverow,c)/circ_std(saverow,c);
        end

        clear noisevalues
        clear circvalues
        clear row

    end %loop c cellline

end %loop a reporter

fig = figure;
fig.Position = [420,285,525,425];

for j = 1:3 %loop j different plots

    marker = 's^dhv<ox>';

    if j < 3 %plot circadian vs noise component, overlay  Bmal1 and Per2

        x = circ(j,:);
        y = noise(j,:);
        x_err = circ_std(j,:);
        y_err = noise_std(j,:);
        ylabeltext = 'Normalized noise (log10)';
        xlabeltext = 'Normalized circadianicity';
        filename = append('mra_scatterplot_circadianicity_vs_noise');

    elseif j == 3 %plot circadian component Bmal1 vs Per2

        fig = figure;
        fig.Position = [420,285,525,425];

        x = circ(1,1:7);
        y = circ(2,1:7);
        x_err = circ_std(1,1:7);
        y_err = circ_std(2,1:7);
        ylabeltext = 'Per2 circadianicity';
        xlabeltext = 'Bmal1 circadianicity';
        filename = append('mra_scatterplot_circadianicity_bmal1_vs_per2');

    end

    %plot errorbars first so that they appear behind the scatterplot markers
    eb(1) = errorbar(x,y,x_err, 'horizontal', 'LineStyle', 'none','HandleVisibility','off'); hold all
    eb(2) = errorbar(x,y,y_err, 'vertical', 'LineStyle', 'none','HandleVisibility','off');
    set(eb, 'LineWidth', 0.75, 'Color','black')

    if j == 1
        s = gscatter(x,y,celllines,color,marker);
        set(s,'MarkerSize',13,{'MarkerFaceColor'},get(s,'Color'));
        legend('Location','northeastoutside','FontSize',16);
        hold on
    elseif j == 2
        z = gscatter(x(:,1:7),y(:,1:7),celllines(1:7,:),color(1:7,:),marker,13,'HandleVisibility','off');
        set(z,'MarkerSize',13,'LineWidth',2,'MarkerFaceColor','w','HandleVisibility','off');
    elseif j == 3
        q = gscatter(x,y,celllines(1:7,:),color,marker,13);
        set(q,'MarkerSize',13,{'MarkerFaceColor'},get(q,'Color'));

        %linear regression

        % Removing rows with missing values in either x or y
        validRows = ~isnan(allcirc_bmal) & ~isnan(allcirc_per);
        regr_x = allcirc_bmal(validRows);
        regr_y = allcirc_per(validRows);

        mdl = fitlm(regr_x,regr_y);
        plot(mdl,'color','none','HandleVisibility','off','LineWidth',1,'Marker','none');
        str =['R^2 = ',sprintf('%.2f',mdl.Rsquared.Ordinary)];
        annotation('textbox',[.15 0.9 0 0],'string',str,'FitBoxToText','on','EdgeColor','black');

        hold off
    end

    clear x
    clear y
    clear x_err
    clear y_err

    ax = gca;

    %specify remaining settings and appearance of the plot
    ylabel(ylabeltext,'FontSize',16,'FontName',font);
    xlabel(xlabeltext,'FontSize',16,'FontName',font);

    box on;
    grid on;
    set(ax,'LineWidth',0.9,'FontName',font,'FontSize',18,'XMinorGrid','off','YMinorGrid','off');

    if j == 2
        %add cutoff line on log-lin scale
        X=(0:0.0001:1);
        X2 = flip(X);
        plot(X,X2,'Handlevisibility','off','Color',[0.3 0.3 0.3],'LineWidth',1);
        hold off
        set(gca,'YScale','log')
        ylim([0.00007 1.3]);
        yticklabels({'0.0001','0.001','0.01','0.1','1'});
        xlim([-0.03 1.03]);
        xticks([0:0.2:1]);
        %save figure
        saveas(fig, [filename, '.svg']);
    elseif j == 3
        xlim([-0.03 1.03]);
        ylim([-0.03 1.03]);
        %save figure
        saveas(fig, [filename, '.svg']);
    end

end %loop j plots

end %function