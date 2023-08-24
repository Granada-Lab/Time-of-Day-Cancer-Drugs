function boxplot_circadian_values

%Carolin Ector, 22.08.2023
%rank circadian values (autocorrelation peak and lag, ridgelength) per cell line

metric = {'lag';'peak';'ridgelength'};
celllines = {'MCF10A';'MCF7';'HCC1806';'MDAMB468';'GIMEN';'SY5Y';'U2OS';'U2OS sKO';'U2OS dKO'};

for m = 1:3 %loop metric

    fig = figure;
    fig.Position = [420,285,525,425];

    %load data
    if m < 3
        path = 'autocorrelation_results.xlsx';
    elseif m == 3
        path = 'cwt_ridgelengths_threshold_halfmax.xlsx';
    end
    sheet_bmal = append('Bmal1_',metric{m});
    sheet_per = append('Per2_',metric{m});

    [bmal] = xlsread(path,sheet_bmal);
    [per] = xlsread(path,sheet_per);

    %no Per2 data for U2OS-KO cell lines
    per(:,8) = NaN;
    per(:,9) = NaN;

    %merge Bmal1 and Per2 data for all cell lines
    data = [bmal;per];

    if m == 1 %exclude SY5Y and U2OS Cry1/2-dKO from ranking due to periods above the circadian range
        data(:,9) = [];
        data(:,6) = [];
        celllines(9) = [];
        celllines(6) = [];
    end

    %sort cell lines by their median lag or peak value (ascending)
    med = median(data,'omitnan');
    [~, sortIdx] = sort(med,'ascend');
    dataAscend = data(:,sortIdx);
    celllinesAscend = celllines(sortIdx,:);

    %plot boxplot
    boxplot(dataAscend,'Labels',celllinesAscend); hold on

    %add datapoints as scatter
    datasize = size(dataAscend);
    x = repmat(1:datasize(:,2),datasize(:,1),1);
    s = scatter(x,dataAscend,36,'k','filled','MarkerFaceAlpha',0.8','jitter','on','jitterAmount',0.15);
    set(s,'MarkerEdgeColor',[0 0 0]);

    if m == 1
        ylabeltext = 'Lag 2^{nd} peak (h)';
    else
        ylabeltext = 'Autocorrelation 2^{nd} peak';
    end

    ylabel(ylabeltext)

    ax = gca;
    set(findobj(gca,'type','line'),'linew',1);
    box on;
    grid on
    set(ax,'YLimitMethod', 'padded','LineWidth',0.9,'FontName','Helvetica Neue','FontSize',18,'XMinorGrid','off','YMinorGrid','off');

end %loop m metric
end %function