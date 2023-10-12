function ToD_mcf10a_overlay_benefit_extraction(date_all,cellline,drug,mcf10_col,color_all,outputorder)

%Carolin Ector, 23.08.2023

%%Function ...:
%%reads Relative ToD Response data (ToD profile) from excel sheet (calculated with function: ToD_growth_plots_parameter_extraction)
%%overlays ToD profiles of mcf10a cell (healthy) and a cancer cell line
%%identifies effect size and times of maximum and minimum benefit:
% - maximum benefit (highest response difference between cancer cell line and mcf10a with higher toxicity against cancer cell line as opposed to mcf10a)
% - minimum benefit (highest response difference between cancer cell line and mcf10a with higher toxicity against mcf10a as opposed to cancer cell line)
% - minimum difference (minimum response difference between cancer cellline and mcf10a)

%Time-of-Day-Cancer-Drugs Manuscript Fig.4h-k

%input: stored in "ToD_overlay_workspace.mat"
% date_all: dates of all ToD experiments
% cellline: names of the cell lines being analysed
% drug: names of drugs used for time-of-day treatments
% color_all: color for the different cell line models
% mcf10_col: column order to load mcf10a data in the drug (d) loop
% outputorder: column order to save results of the different cell lines in the excel sheet

%Define remaining variables
channel = {'Conf';'CellNr'}; %confluency = brightfield channel, cell number = red fluorescent channel from live-imaging
yaxisnames = {'Confluency (%)';'Cell Number'};
timepoints1 = [0;4;8;16;20;24]; %timepoints of treatment for all drugs except cisplatin
timepoints2 = [0;3;6;9;15;18;21;24]; %timepoints of treatment for cisplatin

cc = numel(cellline);

%% Load and normalize data

for a = 1:2 %loop a channel

    outputfile = 'ToD_values_splinesmooth_benefit.xlsx';

    sheet1 = append('responses_',channel{a});
    sheet2 = append('stdresponses_',channel{a});

    excelfile_mcf10a = append('ToD_experiments_results/20221107_ToD_Results_4d_MCF10A.xlsx');
    [data_mcf10] = readmatrix(excelfile_mcf10a,'Sheet',sheet1);
    [std_mcf10] = readmatrix(excelfile_mcf10a,'Sheet',sheet2);

    for c = 1:cc %loop c cell lines

        row = outputorder(:,c);

        if c == 2
            dd = 7;
        elseif c == 10
            dd = 9;
        else
            dd = 8;
        end

        for d = 1:dd %loop d drug

            if c == 2 && d == 8

            else

                if c < 4 && d < 8
                    date = date_all{1};
                elseif c > 3 && c < 7 && d < 8
                    date = date_all{2};
                elseif c > 6 && c < 10 && d < 8
                    date = date_all{3};
                elseif c == 10
                    date = date_all{4};
                end

                if c ~= 10 && c ~= 2 && d == 8
                    numbersToCheck = [1,6,8];
                    if ismember(c, numbersToCheck)
                        date = date_all{6};
                    else
                        date = date_all{5};
                    end
                end

                % Load excel file where data is stored
                inputfile = append('ToD_experiments_results/',date,'_ToD_Results_4d_',cellline{c},'.xlsx');
                [data] = readmatrix(inputfile,'Sheet',sheet1);
                [stdev] = readmatrix(inputfile,'Sheet',sheet2);

                if c == 10
                    col = mcf10_col(:,d);
                elseif c < 10 && d == 8
                    col = 2;
                else
                    col = d+1;
                end

                figure('Position',[1000,918,672,420]);

                for v = 1:2 %loop mcf10a and cancer cell line

                    if v == 1
                        responses = data_mcf10(:,mcf10_col(:,d));
                        err = std_mcf10(:,mcf10_col(:,d)); %standard dev.
                        color = [0.2,0.2,0.2];
                        col2 = 1;
                        marker = 's';
                    else
                        if d == 3 && c == 9 %SUM140PT Alpelisib
                            responses(1:6,:) = 1;
                            err(1:6,:) = 0;
                        else
                            responses = data(:,col);
                            err = stdev(:,col);
                        end
                        color = cell2mat(color_all(c,:));
                        col2 = c+1;
                        marker = 'o';
                    end

                    if d == 8
                        if v == 1
                            timepoints = timepoints1;
                        elseif v == 2 && c ~= 10
                            timepoints = timepoints2;
                        end
                    else
                        timepoints = timepoints1;
                    end

                    %% spline smoothing

                    %Spline interpolation provides a smooth curve that passes through the given data points.
                    splinefit = fit(timepoints,responses,'smoothingspline','SmoothingParam',0.7);

                    % Evaluate spline fit within the range of timepoints
                    xRange = linspace(min(timepoints), max(timepoints), 1000); % Adjust the number of points as desired
                    yfitinterpolate(:,v) = feval(splinefit, xRange);
                    h1 = plot(xRange, yfitinterpolate(:,v)); hold all

                    if c == 10 && v == 1
                        ToDMR_spline{6,d} = max(yfitinterpolate(:,v))-min(yfitinterpolate(:,v)); %save mcf10a result
                    else
                        ToDMR_spline{row,d} = max(yfitinterpolate(:,v))-min(yfitinterpolate(:,v)); %save result cancer cell line
                    end

                    % fine-tune plots
                    set(h1,'LineWidth',2,'color',color,'HandleVisibility','off');
                    errorbar(timepoints,responses,err,marker,'MarkerSize',15,'LineWidth',2,'MarkerFaceColor',color,'Color',color);

                    % calculate maximum difference in responses between cancer cell line and healthy (MCF10A) cell line
                    if v == 2

                        %% method1: dividing the responses of MCF10A from the responses of the cancer cell line
                        %dividing the responses of the cancer cell line with the responses of MCF10A
                        foldchange = yfitinterpolate(:,1)./yfitinterpolate(:,2);

                        %find maximum benefit and corresponding treatment timepoint (considering timepoints after ToD = 2 hours (=row 83))
                        max_benefit{d,col2} = max(foldchange(83:end,:));
                        idx_max_benefit = find(foldchange == max_benefit{d,col2});
                        time_max_benefit{d,col2} = xRange(idx_max_benefit);

                        %find minimum benefit and corresponding treatment timepoint (considering timepoints after ToD = 2 hours (=row 83))
                        min_benefit{d,col2} = min(foldchange(83:end,:));
                        idx_min_benefit = find(foldchange == min_benefit{d,col2});
                        time_min_benefit{d,col2} = xRange(idx_min_benefit);

                        %find timepoint with minimal difference in treatment responses (~1)
                        % Start from row 83 to the end
                        foldchange_subset = foldchange(83:end, 1);
                        % Find the index of the minimal absolute difference from 1
                        [~, idx_min_diff] = min(abs(foldchange_subset - 1));
                        % Assign the minimal difference value to min_diff
                        min_diff{d, col2} = foldchange_subset(idx_min_diff);
                        % Assign the corresponding timepoint to time_min_diff
                        time_min_diff{d, col2} = xRange(82 + idx_min_diff); % Adjust for the offset from starting at row 83

                        %                         %% method2: substracting the responses of MCF10A from the responses of the cancer cell line
                        %                         diff = yfitinterpolate(:,1)-yfitinterpolate(:,2);
                        %
                        %                         %find maximum benefit and corresponding treatment timepoint (considering timepoints after ToD = 2 hours (=row 83))
                        %                         max_benefit{d,col2} = max(diff(83:end,:));
                        %                         idx_max_benefit = find(diff == max_benefit{d,col2});
                        %                         time_max_benefit{d,col2} = xRange(idx_max_benefit);
                        %
                        %                         %find minimum benefit and corresponding treatment timepoint (considering timepoints after ToD = 2 hours (=row 83))
                        %                         min_benefit{d,col2} = min(diff(83:end,:));
                        %                         idx_min_benefit = find(diff == min_benefit{d,col2});
                        %                         time_min_benefit{d,col2} = xRange(idx_min_benefit);

%                         %find timepoint with minimal difference
%                         % Calculate the absolute differences between each element and 0
%                         abs_differences = abs(foldchange(83:end,:) - 1);
%                         % Find the index of the minimum absolute difference
%                         idx_mindiff = find(abs_differences == min(abs_differences));
%                         % Get the data point closest to 0
%                         min_diff{d,col2} = foldchange(idx_mindiff);
%                         time_min_diff{d,col2} = xRange(idx_mindiff);

                        %mark maximum and minimum benefit times in response curve overlay plot
                        xline(time_max_benefit{d,col2}, 'color', [0.12 0.82 0.37], 'LineWidth',20,'LineStyle','-','Alpha',0.2);
                        xline(time_min_benefit{d,col2}, 'color', [1 0 0], 'LineWidth',20,'LineStyle','-','Alpha',0.2);

                    end

                    clear xData
                    clear yData
                    clear foldchange

                end

                clear yfitinterpolate

                % insert line at y=1 (=control)
                yline(1, 'color', [0.4, 0.4, 0.4], 'LineWidth',2,'LineStyle','--','HandleVisibility','off');

                % add legend
                legendentries = {'MCF10A';cellline{c};'Max Benefit';'Min Benefit'};
                hLeg = legend(legendentries);
                lgdOpt = {'FontSize',16,'Orientation','vertical','Location','northeastoutside', ...
                    'FontWeight','normal','EdgeColor','none','Color','#f5f5f5'};
                set(hLeg,lgdOpt{:});

                titletext = drug{d};

                mainTextOpt = {'FontSize',22,'FontName','Helvetica Neue','FontWeight','normal'};
                axOpt = {'XLimitMethod','padded','YLimitMethod','padded','linewidth',1.5,'YGrid','on','XGrid','off','Box','on','Color','none','FontSize',22};

                xlabel('Time of Day (h)',mainTextOpt{:});
                ylabeltext = append('Relative ',yaxisnames{a});
                ylabel(ylabeltext,mainTextOpt{:});
                title(titletext,mainTextOpt{:});

                ax = gca;
                grid on;

                xticks(0:4:24); xlim([-1 25]);
                set(ax,axOpt{:});

                hold off

                %save figure
                filename = append('ToD_plots/',channel{a},'/overlay_cancercellline_mcf10a/',date,'_ToD_',cellline{c},'_',drug{d},'_',channel{a},'_splinesmooth_max_min_benefit');
                saveas(h1, [filename '.svg']);

            end %exclude HCC1806 --> no cisplatin data

        end %loop d drug

        close all

    end %loop c cell lines

    outputsheetnames = {'maxbenefit_','maxbenefittime_','minbenefit_','minbenefittime_','mindifference_','mindifferencetime_','ToDMR_spline_'};
    inputvariables = {max_benefit;time_max_benefit;min_benefit;time_min_benefit;min_diff;time_min_diff;ToDMR_spline};

    for i = 1:numel(inputvariables)

        input = inputvariables{i};

        outputsheet = append(outputsheetnames{i},channel{a});
        if i == 7
            rownames = cell2table({'HCC1143';'HCC1937';'MDAMB468';'HCC1806';'MDAMB231';'MCF10A';'CAL51';'HCC38';'HCC1937_2';'MDAMB436';'SUM149PT'});
            varname = drug;
            input{11,3} = NaN; %no data for SUM149PT Alpelisib
        else
            varname = [cellline(1:1-1); {'MCF10A'}; cellline(1:end)];
            rownames = cell2table(drug);
            input{3,10} = NaN; %no data for SUM149PT Alpelisib
        end
        table = cell2table(input,"VariableNames",varname);
        finaltable = [rownames,table];
        writetable(finaltable,outputfile,'sheet',outputsheet);

        clear input
        clear outputsheet

    end %loop i saving of values

end %loop a channel
end %function
