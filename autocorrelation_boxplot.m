function autocorrelation_boxplot


path = 'autocorrelation_results_cutoff_48h_v2.xlsx';
sheetnames = {"BMAL_period";"PER_period";"BMAL_RI";"PER_RI"};
metric = {"Period";"Power"};

mm = numel(metric);

for m = 1:mm %loop metric

    celllines = {'MCF10A';'MCF7';'HCC1806';'MDAMB468';'GIMEN';'SY5Y';'U2OS';'U2OS CRY1-sKO';'U2OS CRY1/2-dKO'};
    celllines = transpose(celllines);
    color_string = {'[0,0,0]';'[1,0.07,0.65]';'[0.25,0.14,0.91]';'[0.24,0.61,0.98]';'[0.04,0.75,0.11]';'[0.12,0.98,0.12]';'[0.98,0.64,0.05]';'[1.00,0.00,0.00]';'[0.60,0.00,0.00]'};

    %for c = 1:cc %loop cellline

    fig = figure;
    fig.Position = [420,285,525,425];

    if m == 1
        sheet_bmal = append(sheetnames{1});
        sheet_per = append(sheetnames{2});
        %colororder = [3,6,1,7,5,2,4];
    else
        sheet_bmal = append(sheetnames{3});
        sheet_per = append(sheetnames{4});
        %colororder = [3,8,1,6,5,4,8,7,2];
    end

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

    per(:,8) = NaN; %no PER2 data for KO cell lines
    per(:,9) = NaN;

    data = [bmal;per];

    if m == 1
        data(:,9) = [];
        data(:,6) = [];
        celllines(9) = [];
        celllines(6) = [];
    end 

    med = median(data,'omitnan');
    %stdev = std(data,[],1,'omitnan');
    [~, sortIdx] = sort(med,'ascend');
    dataAscend = data(:,sortIdx);
    celllinesAscend = celllines(:,sortIdx);
    color_string2 = transpose(color_string);
    colorAscend = color_string2(:,sortIdx);

    boxplot(dataAscend,'Labels',celllinesAscend); hold on
    
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        %colorrow = colororder(:,j);
        color(j,:) = str2num(colorAscend{:,j});
        patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.3);
    end

    datasize = size(dataAscend);
    x = repmat(1:datasize(:,2),datasize(:,1),1);
    s = scatter(x,dataAscend,36,color,'filled','MarkerFaceAlpha',0.8','jitter','on','jitterAmount',0.15);
    set(s,'MarkerEdgeColor',[0 0 0]);

    if m == 1
        ylabeltext = 'Autocorrelation Period [h]';
    else
        ylabeltext = 'Rhythmicity Index';
    end

    ylabel(ylabeltext)

    ax = gca;
    set(findobj(gca,'type','line'),'linew',1);
    box on;
    grid on
    set(ax,'YLimitMethod', 'padded','LineWidth',0.9,'FontName','Helvetica Neue','FontSize',18,'XMinorGrid','off','YMinorGrid','off');

    %titletext = append(metric{m});
    %title(titletext,'FontSize',18);

end
end

%save figure
%         filetext = append('phase_difference_output_',cellline{c},'_exp1_vs_exp2','_fig4_boxplot.png');
%         folderpath = strcat('Phase_Difference_Plots/');
%         saveas(fig, [folderpath, filetext]);


%     ii = numel(cellline);
%     for i = 1:ii
%         color(:,i) = str2num(colorAscend{:,i});
%     end

%     for j=1:cc
%         boxchart(celllinesAscend(:,j),'categorical',data(:,j),'BoxFaceColor', color(j,:))
%     end

%     for j=1:6
%         color = str2num(colorAscend{j,:});
%         boxplot(dataAscend(:,j),'Labels',celllinesAscend(:,j),'Colors',color); hold all
%     end


%     data = {[bmal;per]}; %boxchart
%     celllines_string2 = transpose(celllines_string);

%boxchart(categorical(celllines_string2),data);

%     h = boxplotGroup(data,'PrimaryLabels', celllines,'Colors','b','GroupType','betweenGroups', ...
%         'GroupLabelType', 'Vertical'); hold all