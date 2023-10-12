function DR_analysis_pipeline(date,cellline,drug,doses,v,colorarray,nr,dosetextexcel)

%Carolin Ector, 11.09.2023

%Time-of-Day-Cancer-Drugs Manuscript Fig.3 + Fig.S3

%Function processes growth data of dose-response (DR) experiments for multiple drugs and cell lines
%runs the following scripts: 
    % DR_growth_plots: Plot growth curves from DR experiments
    % DR_exponential_fit: Fit an exponential function to growth curves from DR experiments
    % DR_response_curve_fit: Compute growth rate inhibition parameters
    % DR_EC50_extraction: Compute tradiational EC50 values

%input: stored in '[date]_DR_workspace.mat'
% date: date of the experiment being analyzed
% cellline: names of the cell lines being analysed
% drug: names of drugs used for drug treatments of the dose-response experiments
% doses: doses of each drug administered to the cells
% v: column order in excel sheet to load data as follows:
    % 1. solvent control
    % 2. lowest dose (dose(1)),
    % 3. intermediate doses (dose(2),...,dose(n-1))
    % 4. highest dose (dose(n))
% colorarray: colors for different drug doses

% Define remaining variables
experiment = str2num(date);
channel = {'CellNr';'Conf'};
dd = numel(drug); %loop drugs
ee = (numel(doses{1})+1); %loop doses, +1 = control
mkdir('DR_plots'); %create folder to save the different figures in
mkdir('DR_results'); %create folder to save the result files in

%% Load and normalize data
for a = 1:2 %loop a channel

    k_all = cell(1,dd);
    CIlow_all = cell(1,dd);
    CIup_all = cell(1,dd);
    Rsq_all = cell(1,dd);

    cc =  numel(cellline);
    if experiment == 2021 && a == 1
        cc = cc-1; % no Cell Number data for MDAMB436
    end

    for c = 1:cc %loop c cell lines

        % Load excel file where data is stored
        path = append('dose_response_experiments_raw_data/',date,'_DR_',cellline{c},'.xlsx');
        sheet1 = append(channel{a},'_1');
        if experiment == 2021
            sheet2 = append(channel{a},'_1'); %no plate duplicates
        else
            sheet2 = append(channel{a},'_2');
        end

        [num1] = xlsread(path,sheet1);
        [num2] = xlsread(path,sheet2);

        for d = 1:dd %loop d drug

            %define order how to load columns (data organization is partially different for the different experiments)
            %goal: load columns serially by increasing doses starting with the control
            if experiment == 20220204
                yvalues = [v(d),v(d)+5,v(d)+4,v(d)+3,v(d)+2,v(d)+1];
            elseif experiment == 20220317
                yvalues = flip([v(d):7:45,3]);
            elseif experiment == 20221026
                %control:
                if d == 9
                    ctrl_1 = 61; ctrl_2 = 61; %plate1; plate2
                elseif d == 7 || d == 8
                    ctrl_1 = 46; ctrl_2 = 25; %plate1; plate2
                elseif d == 4 || d == 5 || d == 6
                    ctrl_1 = 3; ctrl_2 = 25; %plate1; plate2
                else
                    ctrl_1 = 3; ctrl_2 = 3; %plate1; plate2
                end
                %treated:
                w = [4,11,18,26,33,40,47,54,62]; %plate2
                yvalues_1 = [ctrl_1,v(d),v(d)+1,v(d)+2,v(d)+3,v(d)+4,v(d)+5,v(d)+6]; %plate1
                yvalues_2 = [ctrl_2,w(d),w(d)+1,w(d)+2,w(d)+3,w(d)+4,w(d)+5,w(d)+6]; %plate2
            elseif experiment == 2021
                if c == 1 || c == 10
                    yvalues = [3,12,10,8,6,4,13,11,9,7,5];
                else
                    yvalues = [3,14,12,10,8,6,4,13,11,9,7];
                end
            end

            %load x-values (elapsed recording time)
            xdata1=num1(:,2);
            xdata2=num2(:,2);
            meanx = mean([xdata1,xdata2], 2);

            for e = 1:ee %loop e concentrations

                %load and smooth growth data
                if experiment ~= 2021
                    if experiment == 20221026
                        ysmooth1=smooth(num1(:,yvalues_1(e)),0.35,'rloess');
                        ysmooth2=smooth(num2(:,yvalues_2(e)),0.35,'rloess');
                    else
                        ysmooth1=smooth(num1(:,yvalues(e)),0.35,'rloess');
                        ysmooth2=smooth(num2(:,yvalues(e)),0.35,'rloess');
                    end
                    ydata_smooth(:,e) = mean([ysmooth1,ysmooth2], 2);
                    stdev_smooth(:,e) = std([ysmooth1,ysmooth2],[],2);
                else
                    ysmooth1=smooth(num1(:,yvalues(e)),0.35,'rloess');
                    ysmooth2=smooth(num1(:,yvalues(e)),0.35,'rloess');
                    ydata_smooth(:,e) = ysmooth1;
                    stdev_smooth(:,e) = smooth(num1(:,(yvalues(e)+12)),0.35,'rloess');
                end

                % 1. normalize ydata to time of treatment (= row number)
                timeoftreat = 1;
                if experiment == 2021 && c == 1
                    timeoftreat = 20;
                end
                imaginginterval = xdata1(2,1);

                % adjust lengths of datasets to range from time of treatment to 96h post treatment
                row_end = timeoftreat + 96/imaginginterval; % find corresponding row number where time after treatment = 96h
                ysmooth1_2 = ysmooth1(timeoftreat:row_end,:);
                ysmooth2_2 = ysmooth2(timeoftreat:row_end,:);

                y_normx0_1(:,e) = ysmooth1_2./ysmooth1_2(1,:);
                y_normx0_2(:,e) = ysmooth2_2./ysmooth2_2(1,:);
                ydata_normx0(:,e) = mean([y_normx0_1(:,e),y_normx0_2(:,e)], 2);
                if experiment == 2021
                    stdev_normx0(:,e) = stdev_smooth(timeoftreat:row_end,e)./ysmooth1_2;
                else
                    stdev_normx0(:,e) = std([y_normx0_1(:,e),y_normx0_2(:,e)],[],2);
                end

                % save values in cell array for exponential fit
                %                 xData_all{d} = meanx;
                yNorm_all{e,d} = ydata_normx0(:,e);
                yStd_all{e,d} = stdev_normx0(:,e);
                ymin_all{e,d} = min(ydata_normx0(:,e));
                ymax_all{e,d} = median(ydata_normx0(end-3:end,e));

                % 2. normalize ydata_normx0 to control over time
                yval_normctrl_1 = y_normx0_1(:,e)./y_normx0_1(:,1);
                yval_normctrl_2 = y_normx0_2(:,e)./y_normx0_2(:,1);
                ydata_normctrl(:,e) = mean([yval_normctrl_1,yval_normctrl_2], 2);
                if experiment == 2021
                    stdev_normctrl(:,e) = stdev_normx0(:,e)./y_normx0_1(:,1);
                else
                    stdev_normctrl(:,e) = std([yval_normctrl_1,yval_normctrl_2],[],2);
                end

                valuestoclear1 = {'ysmooth1';'ysmooth2'};
                clear(valuestoclear1{:});

                % 3. extract final responses (for sigmoidal curve fit, traditional EC50 extraction)
                finalresponse{e,d} = median(ydata_normctrl(end-3:end,e));
                sd_finalresponse{e,d} = median(stdev_normctrl(end-3:end,e));

            end %loop concentration

            %%plot growth curves (smoothed raw growth, normalized growth to t=0 and relative growth to control)
            DR_growth_plots(a,c,d,date,channel,cellline,drug,doses,colorarray,imaginginterval,meanx,ydata_smooth,stdev_smooth,ydata_normx0,stdev_normx0,ydata_normctrl,stdev_normctrl)

        end %loop d drug

        %exponential fit of growth curves
        DR_exponential_fit(a,c,date,channel,cellline,drug,doses,nr,dosetextexcel,imaginginterval,yNorm_all,yStd_all,ymin_all,ymax_all)

        %sigmoidal curve fit - Hafner et al. 2017 approach
        [parameters, dose_GRinf] = DR_response_curve_fit(a,c,date,experiment,cellline,drug,doses);

        %sigmoidal curve fit - traditional EC50 value extraction
        [parameters] = DR_EC50_extraction(a,c,date,cellline,drug,doses,finalresponse,sd_finalresponse,parameters);

        %save smallest dose elicing a GRinf response for each drug
        growth_GRinf(:,1) = (0:imaginginterval:96);
        growth_GRinf_std(:,1) = (0:imaginginterval:96);
        columnheader{1} = 'Time';
        for d2 = 1:dd
            row = dose_GRinf(2,d2)+1; %+1 since first row in yNorm_all = control 
            growth_GRinf(:,d2+1) = cell2mat(yNorm_all(row,d2));
            growth_GRinf_std(:,d2+1) = cell2mat(yStd_all(row,d2));
            columnheader{1,d2+1} = append(drug{d2},'_',num2str(dose_GRinf(1,d2)));
        end

        close all

        % %create excel file where drug sensitivity parameters will be stored
        valuestosave = {parameters;dose_GRinf;growth_GRinf;growth_GRinf_std};
        outputsheets ={'parameters_';'dose_GRinf_';'growth_GRinf';'growth_GRinf_std'};

        for u = 1:numel(valuestosave)
            if u > 2
                outputfile=append('DR_results/',date,'_Growth_GRinf_',cellline{c},'.xlsx');
            else
                outputfile=append('DR_results/',date,'_DR_Parameters_',cellline{c},'.xlsx');
            end
            outputsheet = append(outputsheets{u},'_',channel{a});
            values = valuestosave{u};
            if u == 1
                output = array2table(values,"VariableNames",drug);
                valuenames = {'GEC50';'GR50';'GRinf';'Hill';'GRAOC';'EC50'};
                rownames = cell2table(valuenames);
                outputtable = [rownames,output];
            elseif u == 2
                outputtable = array2table(values,"VariableNames",drug);
            elseif u == 3 || u == 4
                outputtable = array2table(values,"VariableNames",columnheader);
            end
            writetable(outputtable,outputfile,'sheet',outputsheet);
        end

        valuestoclear2 = {'xdata1';'xdata2';'meanx';'yvalues_1';'yvalues_2';'ydata_smooth';'stdev_smooth';'ydata_normx0';'stdev_normx0';
            'y_normx0_1';'y_normx0_2';'yval_normctrl_1';'yval_normctrl_2';'ydata_normctrl';'stdev_normctrl';
            'finalresponse';'sd_finalresponse';'parameters';'growth_GRinf';'growth_GRinf_std';'columnheader'};
        clear(valuestoclear2{:});

    end %loop c cell lines
end %loop a channel
end %function