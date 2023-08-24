function overlay_detrended_signal_mean(file,celllines,col_bmal,col_per)

%Carolin Ector, 23.08.2023
%overlay detrended signals from the Bmal1 and Per2 luciferase reporters per cell line

%input: stored in overlay_detrended_signal_mean_workspace.mat
% file: excel sheet where detrended circadian time-series data is stored
% celllines: names of the cell lines being analysed
% col_bmal / col_per: columns in the excel sheet where respective time series from the Bmal1-reporter or Per2-reporter is stored for replicate1

reporter = {'Bmal1';'Per2'};
ii = numel(reporter);
cc = numel(celllines);

[num] = xlsread(file,'detrended');
time = num(:,1);

for c = 1:cc %loop c celllines

    %load data for bmal1 and per2 reporter
    rep = 3; %take the first 3 replicates to overlay signals from a single experiment

    if c > 7 %KO celllines have no data of the Per2 reporter
        ii = 1;
    end

    for i = 1:ii %loop i reporter

        for r = 1:rep
            a = 0:(rep-1);
            if i == 1
                Bmal1(:,r) = num(:,col_bmal(c)+a(r));
            else
                Per2(:,r) = num(:,col_per(c)+a(r));
            end
        end
    end

    %% process and plot data

    Bmal1_mean = (mean(Bmal1,2));
    Bmal1_std = (std(Bmal1,[],2));
    if c < 8
        Per2_mean = (mean(Per2,2));
        Per2_std = (std(Per2,[],2));
    end

    x = transpose(time);
    y1 = transpose(Bmal1_mean);
    e1 = transpose(Bmal1_std);

    color1 =  [0.9804,0.6392,0.0510]; %yellow for Bmal1

    if c < 8
        y2 = transpose(Per2_mean);
        e2 = transpose(Per2_std);
        color2 = [0.0627,0.4510,0.6824]; %blue for Per2
    end

    fig = figure;

    hold all;

    if c < 8
        yyaxis left
        patch([x fliplr(x)], [(y1-e1)  (fliplr(y1+e1))], color1, 'FaceAlpha',0.1, 'EdgeAlpha',0, 'HandleVisibility','off');
        plot(x,y1,'LineWidth',2.5,'LineStyle','-','Color',color1);
        ylabel('Detrended Bmal1-Luc signal');
        yyaxis right
        patch([x fliplr(x)], [(y2-e2)  (fliplr(y2+e2))],color2, 'FaceAlpha',0.1, 'EdgeAlpha',0, 'HandleVisibility','off');
        plot(x,y2,'LineWidth',2.5,'LineStyle','-','Color',color2);
        ylabel('Detrended Per2-Luc signal');
    else
        patch([x fliplr(x)], [(y1-e1)  (fliplr(y1+e1))], color1, 'FaceAlpha',0.1, 'EdgeAlpha',0, 'HandleVisibility','off');
        plot(x,y1,'LineWidth',2.5,'LineStyle','-','Color',color1);
        ylabel('Detrended Bmal1-Luc signal');
    end

    hold off

    xlabel('Time (d)','FontSize',20,'FontName','Helvetica Neue');

    titletext = append(celllines{c});
    title(titletext,'FontSize',20,'FontName','Helvetica Neue');

    ax = gca;
    grid on;
    xticks(0:24:144);
    xticklabels({'0','1','2','3','4','5'});

    set(ax,'XLimitMethod','padded','YLimitMethod','padded','linewidth',1.5,'YGrid','on', ...
        'XGrid','off','Box','on','Color','none','FontSize',18,'FontName','Helvetica Neue');

    % add legend
    lgd = legend(reporter);
    set(lgd,'FontSize',15,'Orientation','vertical','Location','northeast','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');

    %save figure
    filetext = append('Detrended_Signal_Overlay_CutOff_48h_',celllines{c},'_mean.svg');
    %     saveas(fig, filetext);

    clear Bmal1

    if c < 8
        clear Per2
    end

end %loop c celllines
end %function
