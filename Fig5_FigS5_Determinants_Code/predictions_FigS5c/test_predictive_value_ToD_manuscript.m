function test_predictive_value_ToD_manuscript(inputfile,datasets,axOpt)

%Carolin Ector, 04.04.2024

%ToD-framework manuscript Fig. S5c

%Function calculates the predictive value of linear regressions shown in ToD manuscript Fig. 5b-e with new data
%loads ToD_MR and predictive values (circadian, growth, drug sensitivity parameters) and splits them into original and new data
%loops through different datasets of predictive value datasets
%loop throigh the different drugs and individual predictive values
%performs linear regression and calculates pearson correlation coefficients for each data-pair
%creates a composite per drug with all the different individual predictive values (eg. all circadian parameters)
%predicts ToD_MR values based on new data
%calculates the mean absolute error between predicted and actual datapoints

%input: stored in "workspace_predictions.mat"
% inputfile: excel file with the different datasets (located on separate excelsheets) and the ToD_MR values
% datasets: predictive values (circadian, growth, drug sensitivity parameters)
% axOpt: settings for figure

test_corrcoefficit = 0;
color_celllines = lines(7);
color_celllines(1,:) = [];

%all linear regressions associated with absolute pearson correlation coefficients >= 0.5 in one figure
fig3 = figure('Position',[1,77,2560,1285]);

for a = 1:numel(datasets) %loop a datasets (circadian / growth / drug sensitivity parameters)

    % Load data from Excel file
    [inputdata1,values] = xlsread(inputfile, datasets{a}); %predictive values
    [inputdata2,drugs] = xlsread(inputfile, 'todmr','A1:H17'); %ToD values

    %process descriptions for rows/columns
    celllines = values(2:end,1);
    values = values(1,2:end);
    drugs = drugs(1,2:end);

    %create empty arrays to store predicted and actual values per cell line (rows) and drug (columns)
    oo = numel(values);
    predictedValues2 = nan(6,7);
    actualValues = nan(6,7);
    x_parameter = nan(6,7);

    %adapt lengths of both data arrays to be equal
    if a >= 3 %no drug sensitivity data for U2OS knockout cell lines
        le = length(inputdata1);
        inputdata2 = inputdata2(1:le,:);
    else
        le = length(inputdata2);
        inputdata1 = inputdata1(1:le,:);
    end

    % Split original and new data
    Original_Data = inputdata1(1:10,:); %Predictors
    New_Data = inputdata1(11:end,:);
    ToDMR = inputdata2(1:10,:); %ToDMR values
    New_ToDMR = inputdata2(11:end,:);

    if a >= 3 %drug sensitivity data only compared with ToD_MR values of the same drugs
        d1 = a-2;
        dd = d1;
    else
        d1 = 1;
        dd = numel(drugs);
    end

    if a < 3
        f=0;
    end

    for d = d1:dd %loop d different drugs

        data1 = nan(6,oo);
        data2 = nan(6,oo);

        cc1 = length(Original_Data);

        predictions_LOOCV = nan(cc1,oo);
        actuals_LOOCV = nan(cc1,oo);

        mae_pred = nan(6,oo);
        mae_loocv = nan(6,oo);

        for o = 1:oo %loop o predictor values

            x = Original_Data(:,o);
            y = ToDMR(:,d);
            x2 = New_Data(:,o);
            y2 = New_ToDMR(:,d);

            % Find rows with any NaN values within each array
            invalidxRows = any(isnan(x),2); % A column vector of logicals for input
            invalidyRows = any(isnan(y),2); % A column vector of logicals for comparison

            % Combine the indicators to identify rows to be marked as NaN in both arrays
            rowsToNaN = invalidxRows | invalidyRows;

            % Mark the identified rows as NaN in both arrays
            % Note: For columns, we use ':' to apply the change to all columns of the identified rows
            x(rowsToNaN,:) = [];
            y(rowsToNaN,:) = [];

            % Find rows with any NaN values within each array
            invalidx2Rows = any(isnan(x2),2); % A column vector of logicals for input
            invalidy2Rows = any(isnan(y2),2); % A column vector of logicals for comparison

            % Combine the indicators to identify rows to be marked as NaN in both arrays
            rowsToNaN2 = invalidx2Rows | invalidy2Rows;

            % Mark the identified rows as NaN in both arrays
            % Note: For columns, we use ':' to apply the change to all columns of the identified rows
            x2(rowsToNaN2,:) = [];
            y2(rowsToNaN2,:) = [];

            if a < 3
                col = d;
            else
                col = 1;
            end

            %% calculate pearson correlation coefficient on original data
            [r(o,col), pValue] = corr(x, y, 'Type', 'Pearson');
            corrcoeff = round(abs(r(o,col)),1);

            % Initialize arrays to hold LOOCV predictions and actual values
            %             predictions_LOOCV = zeros(size(x));
            %             actuals_LOOCV = zeros(size(y));

            for i = 1:length(x) % LOOCV iteration
                % Leave-one-out by excluding the ith observation
                x_train = x([1:i-1, i+1:end], :);
                y_train = y([1:i-1, i+1:end], :);

                x_test = x(i, :);
                y_test = y(i, :);

                % Linear regression on the training set
                mdl1 = fitlm(x_train, y_train);

                if corrcoeff >= 0.5
                    % Predict the outcome for the left-out observation
                    predictions_LOOCV(i,o) = predict(mdl1, x_test);
                    actuals_LOOCV(i,o) = y_test;
                else
                    predictions_LOOCV(i,o) = NaN;
                    actuals_LOOCV(i,o) =  NaN;
                end

            end

            % Calculate the Mean Absolute Error (MAE) for the LOOCV
            mae_LOOCV(o,col) = mean(abs(predictions_LOOCV(:,o) - actuals_LOOCV(:,o)));

            %% Linear regression original data
            % Linear regression
            mdl = fitlm(x, y);

            %% Predictions new data

            % Calculate ToD_MR predictions from new data
            predictedValues = predict(mdl, x2);

            predictedValues = max(predictedValues,0); %no negative ToDMR values possible

            % Calculate MAE per predictive value
            mae(o,col) = mean(abs(predictedValues - y2));

            %save predictedValues and actualValues per cell line
            k = 0;
            n=numel(rowsToNaN2);

            for h = 1:n %loop n celllines
                TF = rowsToNaN2(h,:);
                if TF == 0
                    k = k+1;
                    predictedValues2(h,d) = predictedValues(k,1);
                    actualValues(h,d) = y2(k,1);
                    x_parameter(h,d) = x2(k,1);
                    if corrcoeff >= 0.5
                        data1(h,o) = predictedValues(k,1);
                        data2(h,o) = y2(k,1);
                    else
                        data1(h,o) = NaN;
                        data2(h,o) = NaN;
                    end
                elseif TF == 1
                    predictedValues2(h,d) = NaN;
                    actualValues(h,d) = NaN;
                    x_parameter(h,d) = NaN;
                end
            end

            g = 0;

            if a == 1
                if o ==4%exclude amplitude
                    g = 1;
                end
            end

            if corrcoeff >= 0.5
                if g == 0 %show linear regressions with correlation values >= 0.5 in one figure (all datasets combined)
                    test_corrcoefficit = test_corrcoefficit+1;
                    figure(fig3);
                    subplot(4,5,test_corrcoefficit)

                    hold all

                    % Plot linear regression
                    plot(mdl);

                    % Plot original datapoints
                    scatter(x, y,100,'o','filled','k');

                    %set graph appearance options
                    str = ['r=', sprintf('%.2f', corrcoeff), ', MAE=', sprintf('%.2f', mae(o,col))];
                    if a < 3
                        ylabeltext = append('ToD_MR ',drugs{d});
                    else
                        ylabeltext = append('ToD_MR ',drugs{a-2});
                    end
                    ylabel(ylabeltext,'FontSize',15,'FontName','Helvetica Neue');
                    title(str,'FontSize',15,'FontName','Helvetica Neue')

                    % Plot predicted and actual new data
                    for v = 1:n
                        TF2 = rowsToNaN2(v,:);
                        if TF2 == 0
                            scatter(x_parameter(v,d), actualValues(v,d),80,color_celllines(v,:),'v');
                        end
                    end

                    xlabel(values{:,o},'FontSize',15,'FontName','Helvetica Neue','Interpreter','none');
                    ylabel(ylabeltext,'FontSize',15,'FontName','Helvetica Neue');
                    title(str,'FontSize',15,'FontName','Helvetica Neue');
                    legend off

                    ax = gca;
                    grid on;
                    set(ax,axOpt{:});

                    hold off

                end
            end

            vars2clear = {'x';'y';'x2';'y2';'rowsToNaN';'rowsToNaN2';'invalidxRows';'invalidx2Rows';'invalidyRows';'invalidy2Rows';'predictedValues'};
            clear(vars2clear{:});

        end %loop o preditive values

        mae_pred = mae;
        mae_loocv = mae_LOOCV;

        ooo = oo;

        if a == 1 %exclude amplitude
            data1(:,4)=[];
            data2(:,4)=[];
            mae_pred(4,:)=[];
            mae_loocv(4,:)=[];
            predictedValues2(:,4)=[];
            actualValues(:,4) = [];
            predictions_LOOCV(:,4) = [];
            actuals_LOOCV(:,4) = [];
            ooo = ooo-1;
        end

        %new data
        if a < 3
            Data1_all{d} = {data1};
            Data2_all{d} = {data2};
        else
            Data1_all{a-2} = {data1};
            Data2_all{a-2} = {data2};
        end

        if a >= 3
            s=1;
        else
            s=0;
        end

    end %loop d drugs

    if a < 3

        %Initialize 3D matrices
        x = Data1_all;
        y = Data2_all;

        %cellline-by-cellline
        Data1 = zeros(6, ooo, length(x));
        Data2 = zeros(6, ooo, length(y));

        %Populate 3D matrices with data from cell arrays
        for i = 1:length(x)
            Data1(:,:,i) = cell2mat(x{i});
            Data2(:,:,i) = cell2mat(y{i});
        end


        % Define labels, title, and group names for legend
        Labels = {'Predicted', 'Measured ToD_{MR} value'};
        Title = 'Bland-Altman Plot Across Drugs';
        GroupNames = values;
        if a == 1 %exclude amplitude
            GroupNames(:,4) =[]; %exclude amplitude
        end

        % Define a set of symbols for the 11 measurements (columns)
        Symbols = 'osd^v<>phx+*'; % Ensure there are at least as many symbols as measurements

        % Define colors for the 6 drugs
        % MATLAB default color order is 'brgmcky', extend it if you have more than 7 groups
        Colors = 'brgmck'; % Repeat colors or add new ones to ensure there are enough % Extend this if you have more than 6 drugs

        % Call the BlandAltman function with your 3D data, labels, title, and group names
        [rpc, fig] = BlandAltman(Data1, Data2, Labels, Title, GroupNames, 'symbols', Symbols, 'colors', Colors);

        ax = gca;
        set(ax,axOpt{:});

        figurename3 = append('bland_altmann_plot_',datasets{a},'_todmr_corrcoeff05_final.svg');
        saveas(fig,figurename3);

        clear Data1
        clear Data2

    end

    % Calculate MAE per cell line
    n_percell = sum(~isnan(predictedValues2),2); %number of available drugs for that cell line
    mae_percell = mean(abs(predictedValues2 - actualValues),2,'omitnan');
    mae_percell_std = std(abs(predictedValues2 - actualValues),[],2,'omitnan');

    newcelllines = celllines(11:(11+n-1),:);

    for m = 1:n %n=numel new cell line data
        labels_newcelllines(m,:) = append(newcelllines(m,:),' n=',num2str(n_percell(m,:)));
    end %for m = 1:n

    nonNaNColumns = sum(~all(isnan(predictedValues2), 1));
    xbar = mae_percell;
    ebar = mae_percell_std;
    labels = labels_newcelllines;

    if a < 3

        if nonNaNColumns > 1

            fig2=figure;
            % Creating a bar plot that compares the mean MAE values per cell line
            bar(xbar, 'FaceColor', 'flat');

            % Adding error bars - note the use of 'k' for black color and '.' to denote the style
            hold on; % Keeps the bar graph visible when adding error bars
            errorbar(1:length(xbar), xbar, ebar, 'k', 'linestyle', 'none');

            % Labeling the x-axis with cell line names
            set(gca, 'xtick', 1:length(xbar), 'xticklabel', labels,'FontSize',15);
            ylabel('Mean Absolute Error','FontSize',15,'FontName','Helvetica Neue');
            title(datasets{a},'FontSize',20,'FontName','Helvetica Neue','Interpreter','none');
            hold off

            figurename2 = append('bar_mean_mae_by_cellline_',datasets{a},'_todmr_final.svg');
            saveas(fig2,figurename2);

        end %if a > 3
    end

    %save drug sensitivity results in same array
    if a >= 3
        mae_sens(:,(a-2)) = mae;
        mae_sens_percell(:,(a-2)) = mae_percell;
        mae_LOOCV_sens(:,(a-2)) = mae_LOOCV;
    end %if a >= 3

    outputexcel = append('mae_predicted_vs_actual_todmr.xlsx');
    t_values = cell2table(values');
    t_celllines = cell2table(celllines(11:end,:));
    outputsheet1 = append(datasets{a},'_per_cellline');
    outputsheet2 = append('LOOCV_',datasets{a});

    if a < 3 %circadian values and growth values
        % Save MAE to an Excel file
        t_mae = array2table(mae,'VariableNames',drugs);
        t_final1 = [t_values,t_mae];
        writetable(t_final1, outputexcel,'Sheet',datasets{a});

        t_mae_percell = array2table(mae_percell);
        t_final2 = [t_celllines,t_mae_percell];
        writetable(t_final2, outputexcel,'Sheet',outputsheet1);

        t_mae_LOOCV = array2table(mae_LOOCV,'VariableNames',drugs);
        t_final3 = [t_values,t_mae_LOOCV];
        writetable(t_final3, outputexcel,'Sheet',outputsheet2);
    end %if a < 3

    varstoclear = {'mae';'t_mae';'mae_percell';'t_mae_percell';'mae_LOOCV';'t_mae_LOOCV';'t_final1';'t_final2';'t_final3';'t_values';'predictedValues2';'actualValues';'x';'y';'inputdata1';'inputdata2'};
    clear(varstoclear{:});

end %loop a datasets

if s == 1

    %Initialize 3D matrices
    x = Data1_all;
    y = Data2_all;

    Data1 = zeros(6, ooo, length(x)); %predictions
    Data2 = zeros(6, ooo, length(y));

    %Populate 3D matrices with data from cell arrays
    for i = 1:length(x)
        Data1(:,:,i) = cell2mat(x{i});
        Data2(:,:,i) = cell2mat(y{i});
    end

    % Define labels, title, and group names for legend
    Labels = {'Predicted ToD_{MR} values', 'Actual ToD_{MR} values'};
    Title = 'Bland-Altman Plot Across Drugs';
    GroupNames = values;

    % Define a set of symbols for the 11 measurements (columns)
    Symbols = 'osd^v<>phx+*'; % Ensure there are at least as many symbols as measurements

    % Define colors for the 6 drugs
    % MATLAB default color order is 'brgmcky', extend it if you have more than 7 groups
    Colors = 'brgmckyrbgmcky'; % Repeat colors or add new ones to ensure there are enough % Extend this if you have more than 6 drugs

    % Call the BlandAltman function with your 3D data, labels, title, and group names
    [rpc, fig] = BlandAltman(Data1, Data2, Labels, Title, GroupNames, 'symbols', Symbols, 'colors', Colors);

    figurename3 = append('bland_altmann_plot_drug_sensitivity_todmr_corrcoeff05_final.svg');
    saveas(fig,figurename3);

    % Calculate MAE per cellline
    n_percell2 = sum(~isnan(mae_sens_percell),2); %number of available drugs for that cell line
    for m = 1:n %n=numel new cell line data
        labels_newcelllines2(m,:) = append(newcelllines(m,:),' n=',num2str(n_percell2(m,:)));
    end %for m = 1:n

    xbar = mean(mae_sens_percell,2,'omitnan'); %("1" = drug )
    ebar = std(mae_sens_percell,[],2,'omitnan'); %("1" = drug )
    labels = labels_newcelllines2;
    figurename4 = append('bar_mean_mae_by_cellline_drug_sensitivity_todmr_corrcoeff05_final.svg');

    fig4=figure;
    % Creating a bar plot that compares the mean MAE values per cell line
    bar(xbar, 'FaceColor', 'flat');

    % Adding error bars - note the use of 'k' for black color and '.' to denote the style
    hold on; % Keeps the bar graph visible when adding error bars
    errorbar(1:length(xbar), xbar, ebar, 'k', 'linestyle', 'none');

    % Labeling the x-axis with cell line names
    titletext2 = append('drug_sensitivity');
    set(gca, 'xtick', 1:length(xbar), 'xticklabel', labels,'FontSize',15);
    ylabel('Mean Absolute Error','FontSize',15,'FontName','Helvetica Neue');
    title(titletext2,'FontSize',20,'FontName','Helvetica Neue','Interpreter','none');
    hold off

    saveas(fig4,figurename4);

end

t_values = cell2table(values');

if a >= 3 %collected drug sensitivity values
    %Save MAE to an Excel file
    t_mae = array2table(mae_sens,'VariableNames',drugs);
    t_final1 = [t_values,t_mae];
    writetable(t_final1, outputexcel,'Sheet','drug_sensitivity');

    t_mae_percell = array2table(mae_sens_percell);
    t_final2 = [t_celllines,t_mae_percell];
    writetable(t_final2, outputexcel,'Sheet','drug_sensitivity_per_cellline');

    t_mae_LOOCV = array2table(mae_LOOCV_sens,'VariableNames',drugs);
    t_final3 = [t_values,t_mae_LOOCV];
    writetable(t_final3, outputexcel,'Sheet','LOOCV_drug_sensitivity');
end %if a >= 3

% overview figure compiled of all linear regressions where the correlation coefficient if >= 0.5
han2=axes(fig3,'visible','off');
titletext = append('o = original data, v = actual ToD_MR values',newline,' ');
han2.Title.Visible='on'; han.XLabel.Visible='off'; han.YLabel.Visible='off';
title(han2,titletext,'FontSize',20,'Interpreter','none');
figurename5 = append('linear_regr_predictions_todmr_corcoeff05_final');
saveas(fig3,[figurename5,'.svg']);

end %function