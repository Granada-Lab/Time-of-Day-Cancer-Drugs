function autocorrelation_boxplot(celllines)

%Carolin Ector, 22.08.2023
%rank autocorrelation peak and lag values

%input:
% celllines: names of the cell lines being analysed

metric = {'lag';'peak'};

for m = 1:2 %loop metric

    fig = figure;
    fig.Position = [420,285,525,425];

    %load data
    path = 'autocorrelation_results.xlsx';
    sheet_bmal = append('Bmal1_',metric{m});
    sheet_per = append('Per2_',metric{m});

    [bmal] = xlsread(path,sheet_bmal);
    [per] = xlsread(path,sheet_per);

    %     excl_SY5Y_bmal = [2,5,6]; %periods below 16h or above 50h
    %     excl_SY5Y_per = [4,5,6];
    %     for i = 1:3
    %         nr_b = excl_SY5Y_bmal(:,i);
    %         nr_p = excl_SY5Y_per(:,i);
    %         bmal(nr_b,6) = NaN; %SY5Y, 1st peak not correcty detected for replicate 2 -> exclude
    %         per(nr_p,6) = NaN; %SY5Y, 1st peak not correcty detected for replicate 2 -> exclude
    %     end

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