function butterfly_chart_max_min_benefit_by_drug(file,excelsheets,celllines,channel,drugs)

%Carolin Ector, 23.08.2023

%create butterfly charts to evaluate the extend of the difference in responses between the cancer cell line and MCF10A 
%at the respective maximum and minimum ToD treatment times --> categorized by the individual drugs
% - maximum benefit (highest response difference between cancer cell line and mcf10a with higher toxicity against cancer cell line as opposed to mcf10a)
% - minimum benefit (highest response difference between cancer cell line and mcf10a with higher toxicity against mcf10a as opposed to cancer cell line)

%input: stored in ToD_analysis_workspace.mat
% file: excel file where maximum and minumum benefit foldchanges and times are stored
% excelsheets: sheet names for foldchanges in responses at maximum and minimum benefit
% celllines: names of the cell lines being analysed
% channel: either brightfield (confluency) or fluorescent (cell number) channel from live-imaging, in the manuscript we focus on the fluorescent channel 
% drugs: names of drugs the cell lines have been treated with

dd = numel(drugs);
ii = numel(excelsheets);
celllines(end-3,:) = []; %remove HCC1937_2

for a = 1:2 %loop a channel

    for d = 1:dd %loop d drug

        fig = figure;

        for i = 1:ii %loop i input

            %load data
            inputsheet = append(excelsheets{i},channel{a});
            [alldata] = readmatrix(file,'sheet',inputsheet);

            %save data from all celllines per drug
            data(i,:) = rmmissing(1-alldata(d,:));

            if d == 8
                cellline = {'HCC1143';'MDAMB468';'CAL51';'HCC38';'MDAMB231';'MDAMB436';'SUM149PT';'HCC1937'};
            else
                cellline = celllines;
            end
            
            %save data from all celllines per cell drug
            %             data(i,:) = alldata(d,:); %
            %             x = [1:1:numel(cellline)];

        end %loop i input

        data(:,end-3) = []; %remove HCC1937_2
        x = (1:1:numel(cellline));

        %sort y-axis by highest to lowest maximum benefit
        [y2, sortIdx] = sort(data(2,:),'ascend');
        y1 = data(1,:);
        y1 = y1(:,sortIdx);
        cellline = cellline(sortIdx,:);

        % Create the butterfly chart
        barh(x, y2, 'FaceColor', [0.549 0.9569 0.7373], 'BarWidth', 0.6,'LineWidth',0.8);  %First set of bars
        hold on;
        barh(x, y1, 'FaceColor', [0.9843 0.5804 0.5804], 'BarWidth', 0.6,'LineWidth',0.8);  % Second set of bars (negative values)
        hold off;

        % Adjust axis limits and labels
        ylimits = size(data,2)+0.5;
        ylim([0.5, ylimits]);  
        yticks(1:1:numel(cellline));
        yticklabels(cellline);

        xlabeltext = append('Foldchange ToDMR  Cancer vs. MCF10A ',newline,drugs{d});
        xlabel(xlabeltext);

        % Title and legend
        title('Chronotherapeutic index','FontSize',22);
        lgd = legend('Maximum benefit', 'Minimum benefit');
        set(lgd,'Location','southeast');

        ax = gca;
        set(ax,'XLimitMethod','padded','linewidth',1.5,'Color','w','FontSize',22,'FontName','Helvetica Neue');

        %save figure
        filename = append('butterfly_max_min_ToD_effect_',drugs{d},'_',channel{a});
        saveas(fig, [destination filename '.svg']);

        valuestoclear = {'data';'drug';'sortIdx';'x';'y1';'y2'};
        clear(valuestoclear{:});

    end %loop c cellline / loop d drug
end %loop a metric
end %function