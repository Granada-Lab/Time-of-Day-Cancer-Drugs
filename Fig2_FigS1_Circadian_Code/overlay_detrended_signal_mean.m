function overlay_detrended_signal_mean(file,celllines,col_bmal,col_per,rep_bmal,rep_per)

%Carolin Ector, 23.08.2023

%Function overlays detrended signals from the Bmal1 and Per2 luciferase reporters per cell line

%Time-of-Day-Cancer-Drugs Manuscript Fig. S1c

%input: stored in overlay_detrended_signal_mean_workspace.mat
% file: excel sheet where detrended circadian time-series data is stored
% celllines: names of the cell lines being analysed
% col_bmal / col_per: columns in the excel sheet where respective time series from the Bmal1-reporter or Per2-reporter is stored for replicate1

reporter = {'Bmal1';'Per2'};
cc = numel(celllines);

[num] = xlsread(file,'detrended');
time = num(:,1);

for c = 1:cc %loop c celllines

    if c == 8 || c == 9  %KO celllines have no data of the Per2 reporter
        ii = 1;
    else
        ii = 2;
    end

    for i = 1:ii %loop i reporter

        if i == 1
            rep = rep_bmal; 
        else
            rep = rep_per;
        end

        for r = 1:rep
            a = 0:(rep-1);
            if i == 1
                Bmal1(r,:) = num(:,col_bmal(c)+a(r));
            else
                Per2(r,:) = num(:,col_per(c)+a(r));
            end
        end
    end

    %% process and plot data

    Bmal1_mean = mean(Bmal1,1);
    Bmal1_std = std(Bmal1,[],1);
    if c == 8 || c == 9 
    else
        Per2_mean = mean(Per2,1);
        Per2_std = std(Per2,[],1);
    end

    x = transpose(time);
    y1 = Bmal1_mean;
    e1 = Bmal1_std;

    color1 =  [0.9804,0.6392,0.0510]; %yellow for Bmal1

    if c == 8 || c == 9 
    else
        y2 = Per2_mean;
        e2 = Per2_std;
        color2 = [0.0627,0.4510,0.6824]; %blue for Per2
    end

    fig = figure;

    hold all;

    if c == 8 || c == 9 
        patch([x fliplr(x)], [(y1-e1)  (fliplr(y1+e1))], color1, 'FaceAlpha',0.2, 'EdgeAlpha',0, 'HandleVisibility','off');
        plot(x,y1,'LineWidth',2.5,'LineStyle','-','Color',color1);
        ylabel('Detrended Bmal1-Luc signal');
        lgd = legend(reporter{1});
    else
        yyaxis left
        patch([x fliplr(x)], [(y1-e1)  (fliplr(y1+e1))], color1, 'FaceAlpha',0.2, 'EdgeAlpha',0, 'HandleVisibility','off');
        plot(x,y1,'LineWidth',2.5,'LineStyle','-','Color',color1);
        ylabel('Detrended Bmal1-Luc signal');
        leftAxisLim = max(abs([y1-e1, y1+e1]));
        ylim(leftAxisLim * [-1, 1]);
        leftAxis = gca;
        leftAxis.YColor = color1;
        yyaxis right
        patch([x fliplr(x)], [(y2-e2)  (fliplr(y2+e2))],color2, 'FaceAlpha',0.2, 'EdgeAlpha',0, 'HandleVisibility','off');
        plot(x,y2,'LineWidth',2.5,'LineStyle','-','Color',color2);
        ylabel('Detrended Per2-Luc signal');
        lgd = legend(reporter);
        rightAxisLim = max(abs([y2-e2, y2+e2]));
        ylim(rightAxisLim * [-1, 1]);
        rightAxis = gca;
        rightAxis.YColor = color2;
    end

    hold off

    xlabel('Time (d)','FontSize',20,'FontName','Helvetica Neue');

    titletext = append(celllines{c});
    title(titletext,'FontSize',20,'FontName','Helvetica Neue');

    ax = gca;
    grid on;
    xticks(0:24:144);
    xticklabels({'0','1','2','3','4','5'});

    set(ax,'XLimitMethod','padded','linewidth',1.5,'YGrid','on', ...
        'XGrid','off','Box','on','Color','none','FontSize',18,'FontName','Helvetica Neue');

    % add legend
    set(lgd,'FontSize',15,'Orientation','vertical','Location','northeast','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');

    %save figure
    filetext = append('Detrended_Signal_Overlay_CutOff_48h_',celllines{c},'_mean.svg');
    saveas(fig, filetext);

    clear Bmal1
    clear lgd

    if c == 8 || c == 9 
    else
        clear Per2
    end

end %loop c celllines
end %function
