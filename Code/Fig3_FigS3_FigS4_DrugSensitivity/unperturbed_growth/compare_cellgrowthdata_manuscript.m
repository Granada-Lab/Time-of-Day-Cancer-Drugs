function compare_cellgrowthdata_manuscript(cellline,colorall)

metric = {'Conf';'CellNr'};
cc = numel(cellline);

for c = 1:cc

    for a = 1:2

%         inputfile1 = append('Results_Growth_Data_48_well_DR',metric{a},'.xlsx');
% 
%         [initialvalues{a}] = readmatrix(inputfile1,'sheet','initialvalue');
%         [doublingtimes{a}] = readmatrix(inputfile1,'sheet','doublingime');
%         [growthrates{a}] = readmatrix(inputfile1,'sheet','growthrate');

        inputfile2 = append('Growth_Data_48_well_Selection_Manuscript.xlsx');

        sheet = append(cellline{c});
        [num] = readmatrix(inputfile2,'sheet',sheet);

        time1 = num(:,1);
        %take data from the first 4 days to make different experiments with different recording lengths comparable\
        [~,endrow]=min(abs(time1-96)); %find closest value to 96 in time array

        %define time
        num2 = num(1:endrow,:);
        Time=num2(:,1);
        time_all{:,c} = Time;

        ii = size(num,2)-1;
        for i = 1:ii %loop i replicates
            rawdata = num2(:,i+1);
            smootheddata(:,i) = smooth(rawdata,0.35,'rloess');
            normalizeddata(:,i) = smootheddata(:,i)./smootheddata(1,i);
        end

        %average all replicates per cell line
        mean_smootheddata(:,a) = mean(smootheddata,2,'omitnan');
        std_smootheddata(:,a) = std(smootheddata,[],2,'omitnan');
        mean_normalizeddata(:,a) = mean(normalizeddata,2,'omitnan');
        std_normalizeddata(:,a) = std(normalizeddata,[],2,'omitnan');

        if a == 2
            mean_normalizedall{:,c} = mean_normalizeddata(:,a);
            std_normalizedall{:,c} = std_normalizeddata(:,a);
        end

        varstoclear1 = {'smootheddata';'normalizeddata'};
        
    end

    fig1 = figure(1);

    %add standard deviation as shaded area to plot
    x = transpose(Time);
    y1 = transpose(mean_smootheddata(:,2)); %cell number
    error1 = transpose(std_smootheddata(:,2)); %error cell number
    y2 = transpose(mean_smootheddata(:,1)); %confluence
    error2 = transpose(std_smootheddata(:,1)); %error confluence
    
    patch([x, fliplr(x)], [y1 - error1, fliplr(y1 + error1)],'r','FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'HandleVisibility', 'off');

    colororder({'r','b'})

    hold on;
    yyaxis left
    plot(x,y1,'Color','r','LineWidth',2);
    ylabel('Cell Number');
    yyaxis right
%     ylim([])
    plot(x,y2,'Color','b','LineWidth',2);
    patch([x, fliplr(x)], [y2 - error2, fliplr(y2 + error2)],'b','FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'HandleVisibility', 'off');
    ylabel('Confluence');
    yyaxis right
%     ylim([])

    xlabel('Time [hours]');
    title(cellline{c});

    %create legend
%     legend( b, firstlegendentry, mean_parameters, 'Location', 'NorthWest', 'Interpreter', 'tex','color','w');
    hold off

    % Set the remaining axes and box properties
    ax = gca;
    grid on;
    xticks(0:24:120);
    set(ax,'XLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
        'XGrid','off','Box','on','Color','w','FontSize',16);

    %save figure
    filetext1 = append('exponential_fit_plots/48_well_overlay_cellnr_conf/',cellline{c});
%     saveas(fig1, [ filetext1, '.svg']);

    varstoclear = {'x';'y';'error';'mean_smootheddata';'std_smootheddata';'mean_normalizeddata';'std_normalizeddata';'smootheddata';'normalizeddata'};
    clear(varstoclear{:});
    clf(fig1,'reset');
end

fig2 = figure(2);
cellline2 = {'MCF10A';'HCC1143';'HCC1937';'HCC38';'MDAMB468';'HCC1806';'SUM149PT';'CAL51';'MDAMB231';'MDAMB436'};
dd = numel(cellline2);
colorder = [15,1,3,4,5,6,8,10,11,12];

for d = 1:dd

    column = colorder(:,d);

    x = transpose(cell2mat(time_all(:,column)));
    y = transpose(cell2mat(mean_normalizedall(:,column)));
    error1 = transpose(cell2mat(std_normalizedall(:,column)));

    color = str2num(cell2mat(colorall(column,:)));

    hold all;
    patch([x, fliplr(x)], [y - error1, fliplr(y + error1)],color,'EdgeColor', 'none','FaceAlpha', 0.1, 'HandleVisibility', 'off');
    plot(x,y,'Color',color,'LineWidth',3);

end

ylabel('Normalized Cell Number');
% ylabel('Normalized Conflunecy');
xlabel('Time [hours]');

%create legend
legend(cellline2, 'Location', 'NorthWest', 'Interpreter', 'tex','color','w');

% Set the remaining axes and box properties
ax = gca;
grid on;
xticks(0:24:120);
set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
    'XGrid','off','Box','on','Color','w','FontSize',20);

%save figure
filetext2 = append('norm_cell_nr_48_well_overlay_all_celllines');
% saveas(fig2, [ filetext2, '.svg']);

varstoclear = {'x';'y';'error'};
clear(varstoclear{:});
end