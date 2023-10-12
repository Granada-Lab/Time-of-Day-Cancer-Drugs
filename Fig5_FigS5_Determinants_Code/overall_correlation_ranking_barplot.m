function overall_correlation_ranking_barplot

%Carolin Ector, 23.08.2023

% Function ranks overall correlations between each metric and the maximum range in ToD responses (ToDMR)
    % ranks by metric or by drug
    % ranks absolute or relative average correlations

%Time-of-Day-Cancer-Drugs Manuscript Fig. 5g

%input: stored in "overall_linear_correlation_ranking_workspace.mat"
% labels = labels for x-axis
% mean = mean absolute linear correlation values per drug or metric
% stderr = standard error of absolute linear correlation values per drug or metric
% data tables are stored in corresponding Figshare folder

%load data 
[bmal] = readtable('Data_Fig5/Correlation_Circadian_Bmal1_maxrange_spline.csv');
[drug] = readtable('Data_Fig5/Correlation_DrugResponse_maxrange_spline.csv');
[growth] = readtable('Data_Fig5/Correlation_GrowthProperties_maxrange_spline.csv');

%process drug table to make it compatible with other two tables
drug1 = rows2vars(drug); %transpose
variableNames = drug1{1, :}; % Extract the first row as variable names
variableNames{3} = 'x5FU';
drug1(1, :) = []; % Remove the first row from the table
% Specify the text to add to each variable name
textToAdd = '_maxrange_spline'; % Specify the text to add to each variable name
% Create new variable names by adding text
newVariableNames = cellfun(@(x) [x,textToAdd], variableNames, 'UniformOutput', false);
newVariableNames{1} = 'Var1';
% Set the variable names for the table
drug1.Properties.VariableNames = newVariableNames; % Set the variable names for the table

%sort tables in same order based on bmal1
variableOrder = bmal.Properties.VariableNames;
% Rearrange the columns in each table to match the variable order
bmal2 = bmal(:, variableOrder);
drug2 = drug1(:, variableOrder);
growth2 = growth(:, variableOrder);

%convert to array
bmal3 = table2array(bmal2(1:end,2:end));
drug3 = cell2mat(table2array(drug2(1:end,2:end)));
growth3 = table2array(growth2(1:end,2:end));

%merge data
alldata{1} = vertcat(bmal3,drug3,growth3);
alldata{2} = abs(alldata{1});

%extract names of drugs and parameters for x tick labels
labels{1} = variableOrder(:,2:end)';
labels{2} = vertcat(bmal2.Var1,drug2.Var1,growth2.Var1);

method = {'rel';'abs'};
ranking = {'drug';'parameter'};

for j = 1:2 %loop j relative and absolute values

    %calculate the mean and standard error of the data
    meandata{1} = (mean(alldata{j},1))';
    meandata{2} = mean(alldata{j},2);
    stderr{1} = (std(alldata{j},[],1)./size(alldata{j},2))';
    stderr{2} = std(alldata{j},[],2)./size(alldata{j},1);

    for i = 1:2 %loop i values for ranking (either by drug or metric)

        %load data, i=1: by drug, i=2: by parameter
        labelnames = labels{i};
        y = meandata{i};
        stddata = stderr{i};

        figure;

        x = 1:1:(numel(y));

        %sort y-axis by highest to lowest absolute correlation
        [y, sortIdx] = sort(y,'descend');
        err_ranked = (stddata(sortIdx,:));
        %         err_0(1:numel(y),1) = 0;
        label = labelnames(sortIdx,:);

        hold all;
        h = bar(x,y,'FaceColor', 'r');
        errorbar(x,y,err_ranked,err_ranked,'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
        ylabeltext = append('average ',method{j},' correlation');
        ylabel(ylabeltext);
        hold off

        ax = gca;
        box on;
        grid on;
        if j == 2
            ylim([0.2,Inf]);
        end
        set(ax, 'XTick', 1:1:numel(y),'XTickLabel',label);
        set(ax,'LineWidth',1.5,'FontSize',18,'XMinorGrid','off','YMinorGrid','off','YLimitMethod', 'padded');

        filename = append('average_',method{j},'_correlation_ToD_values','_per_',ranking{i},'.svg');
        saveas(h,filename);

        varstoclear1 = {'labelnames';'x';'y';'staddta'};
        clear(varstoclear1{:})

    end %loop i values for ranking

    varstoclear2 = {'meandata';'stderr'};
    clear(varstoclear2{:})

end %loop j relative and absolute values
end %function
