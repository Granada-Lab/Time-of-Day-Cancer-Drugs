function clustergrams_drug_sensitivity(file,values,values_all,color_crosscorr,cellline,drug,clustergram_colors,displayrange)

%Carolin Ector, 25.08.2023

%Function clusters drug sensitivity values of different cell lines and drugs

%Time-of-Day-Cancer-Drugs Manuscript Fig. 3k-n

%input: stored in 'workspace_DR_clustergrams.mat'
% file: excel file where drug sensitivity values are stored
% values: names of the drug sensitivity values
% cellline: names of the cell lines to be clustered
% drugs: names of the different drugs for the sensitivity values
% clustergram_colors: colors for the different clustergrams
% displayrange: color bar limit for the clustergram

%load data
[data] = xlsread(file,'sensitivity_values');
[GRinf] = xlsread(file,'doses_GRinf'); 
[maxdose] = xlsread(file,'max_tested_dose');

%% clustergram cross correlations sensitivity values: Fig.3n

%columns and rows of the different values and cell lines in the excel file
columns_all = (1:9:46); %initial column for respective drug sensitivity value

for a = 1:numel(values_all)

    val = values_all{a};
    col = columns_all(:,a);

    dataCell{a} = data(:,(col:col+8));

end %loop a values_all

% Calculate correlation coefficients and p-values for each pair of matrices
numDatasets = numel(dataCell);
corrcoeffs = zeros(numDatasets);
pValues = zeros(numDatasets);

for i = 1:numDatasets
    for j = 1:numDatasets
        % Skip self-comparisons
        if i ~= j
            % Flatten the matrices to vectors for correlation calculation
            vector1 = dataCell{i}(:);
            vector2 = dataCell{j}(:);

            % Find common non-NaN indices
            validIndices = ~isnan(vector1) & ~isnan(vector2);

            % Check if vectors have any common non-NaN values
            if any(validIndices)
                % Use only common non-NaN values for correlation calculation
                vector1 = vector1(validIndices);
                vector2 = vector2(validIndices);

                % Calculate correlation coefficients and p-values
                [r, p] = corr(vector1, vector2,'type','Pearson');

                % Store results
                corrcoeffs(i, j) = r;
                pValues(i, j) = p;
            else
                % Set correlation to NaN when no common non-NaN values
                corrcoeffs(i, j) = NaN;
                pValues(i, j) = NaN;
            end
        else
            % Set diagonal elements to 1 (self-comparison)
            corrcoeffs(i, j) = 1;
            pValues(i, j) = 0;  % p-value is 0 for self-comparison
        end
    end
end

% Display p-values (optional)
disp('Correlation Coefficients:');
disp(corrcoeffs);

disp('P-Values:');
disp(pValues);

% Display p-values (optional)
disp('Correlation Coefficients:');
disp(corrcoeffs);

disp('P-Values:');
disp(pValues);

cgo = clustergram(corrcoeffs,'colormap',color_crosscorr);
set(cgo,'RowLabels',values_all,'ColumnLabels',values_all,'AnnotColor','k','AnnotPrecision',3);
set(cgo,'ColumnLabelsRotate',90,'RowLabelsRotate',0);

clear dataCell

valuestosave = {corrcoeffs;pValues};
outputsheets ={'correlationcoeffs';'pvalues'};
outputfile = 'cross_correlation_drug_sensitivity_parameters.xlsx';

for u = 1:numel(valuestosave)
    output = array2table(valuestosave{u},"VariableNames",values_all);
    rownames = cell2table(values_all);
    outputtable = [rownames,output];
    writetable(outputtable,outputfile,'sheet',outputsheets{u});
end

%% individual clustergrams per sensitivity value: Fig.3k-m

%columns and rows of the different values and cell lines in the excel file
columns = (10:9:46); %initial column for respective drug sensitivity value
rows = [3,4,6]; %row for each cell line

for b = 1:numel(values)

    val = values{b};
    col = columns(:,b);

    for c = 1:3 %collect data for all three cell lines
        row = rows(:,c);
        allvalues(c,:) = data(row,(col:col+8));
    end

    allvalues(:,[6,8]) = []; %exclude Alpelisib and Adavosertib (no data for MDAMB468)

    %account for different dose ranges of the drugs:
    %normalize GEC50 and GR50 values to closest dose eliciting a GRinf response
    %normalize GRAOC values to highest tested dose
    if b < 4
        for c = 1:3
            if b < 3
                allvalues(c,:) = (allvalues(c,:)./GRinf(c,:))*100;
            elseif b == 3
                allvalues(c,:) = allvalues(c,:)./maxdose;
            end
        end
    end

    Y = transpose(allvalues);
    cgo = clustergram(Y,'Colormap',clustergram_colors{b});
    set(cgo,'RowLabels',drug,'ColumnLabels',cellline,'AnnotColor','k','AnnotPrecision',3);
    set(cgo,'ColumnLabelsRotate',90,'RowLabelsRotate',0,'Displayrange',displayrange(b,:));
    addTitle(cgo,val);

    clear allvalues

end %loop b values
end %function
