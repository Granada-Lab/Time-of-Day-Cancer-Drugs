% Normalize Excel Data and plot ToD Response Curves

function Plot_DR_Exp_Fit_Plate_By_Plate(folder,date,cellline,drug,doses,v,dosetextexcel)

%Fits an exponential function to growth curves, plate by plate
%INPUT
% date = date of experiment
% v = initial column for respective drug (=CTRL)
% drug = drug as string in cell (8x1 cell = 8 drugs)
% cellline = all cellline as string in cell
% doses = dose range as number in cell (format for 8 drugs: 8x1 cell, each cell containing character array of doses)
%OUTPUT
% k-value from fit
% Ymax value from growth curve
% Relative Growth Factor(RGF) = (relative k) * (relative ymax)

%% Define variables
experiment = str2num(date);
metric = {'CellNr';'Conf'};

aa = numel(cellline);
bb = numel(metric);
cc=numel(drug);
dd = (numel(doses{1})+1); %+1 = control

y_normx0_1_all = cell( dd, cc );
y_normx0_2_all = cell( dd, cc );
ymax_median1 = cell( dd, cc );
ymax_median2 = cell( dd, cc );

%% Load and normalize data
for a = 5:5 %loop a cell lines

    % Load excel file where data is stored
    inputfile = append(folder,date,'_DR_',cellline{a},'.xlsx');

    for b = 1:1 %loop b metric

        %create excel file where data will be stored
        outputfile=append(folder,date,'_Results_ExpFit_high_seed_',cellline{a},'_',metric{b},'.xlsx');

        if b==1
            yaxisname = 'Normalized Cell Number';
        else
            yaxisname = 'Normalized Confluency';
        end

        % Load excel file where data is stored

%         sheet1 = append(metric{b},'_1');
%         sheet2 = append(metric{b},'_2');
        sheet1 = append(metric{b},'_3');
        sheet2 = append(metric{b},'_4');

        [num1] = xlsread(inputfile,sheet1);
        [num2] = xlsread(inputfile,sheet2);

        %take data from the first 4 days to make different experiments with
        %different recording lengths comparable
        if experiment == 20220317
            endrow = find(num1(:,2)==97.5);
        elseif experiment == 20221026
            endrow = find(num1(:,2)==98);
        end

        if experiment == 20220317 || experiment == 20221026
            num1 = num1(1:endrow,:);
            num2 = num2(1:endrow,:);
        end

        %if experiment == 20220317 && b == 3 && a == 1
        %num1(1:8,:) = [];
        %num2(1:8,:) = [];
        %end

        fitresult1 = cell( dd, cc ); %this is saved to the workspace in an answer array and it contains the equation and the calculated values
        fitresult2 = cell( dd, cc );
        gof1 = struct( 'sse', cell( dd, cc ), 'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] ); %gof is somehow not saved to the workspace
        gof2 = struct( 'sse', cell( dd, cc ), 'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
        all_k1 = cell(dd, cc );
        all_k2 = cell(dd, cc );
        allymax1 = cell(dd, cc );
        allymax2 = cell(dd, cc );
        relativegrowthrates1 = cell(dd, cc );
        relativegrowthrates2 = cell(dd, cc );
        Rsq1 = cell(dd, cc );
        Rsq2 = cell(dd, cc );

        for c = 1:cc %loop c drug

            %order: CTRL - lowestdose, dose2, ... highestdose)
            if experiment == 20220204
                yvalues = [v(c),v(c)+5,v(c)+4,v(c)+3,v(c)+2,v(c)+1];
            elseif experiment == 20220317
                yvalues = flip([v(c):7:45,3]);
            elseif experiment == 20221026
                %CONTROL
                if c == 9
                    ctrl_1 = 61; ctrl_2 = 61; %plate1; plate2
                elseif c == 7 || c == 8
                    ctrl_1 = 46; ctrl_2 = 25;
                elseif c == 4 || c == 5 || c == 6
                    ctrl_1 = 3; ctrl_2 = 25;
                else
                    ctrl_1 = 3; ctrl_2 = 3;
                end
                %TREATED
                v = [4,11,18,25,32,39,47,54,62]; %plate1
                w = [4,11,18,26,33,40,47,54,62]; %plate2
                yvalues_1 = [ctrl_1,v(c),v(c)+1,v(c)+2,v(c)+3,v(c)+4,v(c)+5,v(c)+6];
                yvalues_2 = [ctrl_2,w(c),w(c)+1,w(c)+2,w(c)+3,w(c)+4,w(c)+5,w(c)+6];
            elseif experiment == 2021
                if c == 1
                    yvalues = [2,11,9,7,5,3,12,10,8,6,4];
                else
                    yvalues = [2,13,11,9,7,5,3,12,10,8,6];
                end
            end

            %load x-values
            if experiment ~= 2021
                t = 2;
            else
                t = 1;
            end

            Time=num1(:,t);

            for i = 1:dd %loop i concentrations > pre-process data

                % load and smooth raw data doses
                if experiment == 20221026
                    ysmooth1=smooth(num1(:,yvalues_1(i)));
                    ysmooth2=smooth(num2(:,yvalues_2(i)));
                else
                    ysmooth1=smooth(num1(:,yvalues(i)));
                    ysmooth2=smooth(num2(:,yvalues(i)));
                end

                % 1. normalize ydata to timepoint 0
                y_normx0_1 = ysmooth1./ysmooth1(1,:);
                y_normx0_2 = ysmooth2./ysmooth2(1,:);

                % save normalized to x=0 y-values
                y_normx0_1_all{i,c} = y_normx0_1;
                y_normx0_2_all{i,c} = y_normx0_2;
                ymin1{i,c} = min(y_normx0_1);
                ymin2{i,c} = min(y_normx0_2);

%                 %exclude datapoints if there is an initial drop in cell
%                 %numbers (eg. within the first 24 hours, to improve fitting
%                 %quality)
% 
%                 for p = 1:2
%                     if p == 1
%                         data = y_normx0_1;
%                     elseif p == 2
%                         data = y_normx0_2;
%                     end
% 
%                     minvalue = min(data(1:18,:)); %min within first 36h
%                     row_min = find(data == minvalue);
%                     data((1:row_min),:) = NaN;
%                     if p == 1
%                         y_normx0_1_all{i,c} = data;
%                     elseif p == 2
%                         y_normx0_2_all{i,c} = data;
%                     end
% 
%                 end

                %median ymax value from last 3 datapoints
                ymax_median1{i,c} = median(y_normx0_1(end-3:end));
                ymax_median2{i,c} = median(y_normx0_2(end-3:end));

            end %loop concentration

            %% Fit in a loop

            for j = 1:2 % loop well 1 and well

                fig = figure;
% 
                fig.Position = [500,64,940,733];
%                 fig.InnerPosition = [500,64,940,733];
%                 fig.OuterPosition = [500,64,940,812];

                for d = 1:dd

                    if j == 1
                        yMax=ymax_median1{d,c};
                        yMin = min(cell2mat(ymin1(:,c)));
                        yNorm = y_normx0_1_all{d,c};
%                         Plate = 'Plate 1';
                        Plate = sheet1;
                    else
                        yMax=ymax_median2{d,c};
                        yMin = min(cell2mat(ymin2(:,c)));
                        yNorm = y_normx0_2_all{d,c};
%                         Plate = 'Plate 2';
                        Plate = sheet2;
                    end

                    %define y-axis ranges
                    yRange_max = ymax_median1{1,c} + (ymax_median1{1,c}*0.15);
                    if yMin < 1
                        yRange_min = yMin - (yMin*0.15);
                    else
                        yRange_min = 0.5;
                    end

                    [xData, yData] = prepareCurveData(Time, yNorm);

                    % Set up fittype and options.
                    ft = fittype( 'a*exp(x*k)', 'independent', 'x', 'dependent', 'y' );
                    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                    opts.Lower = [0.5 -0.1]; % [Ymin k]
                    opts.Robust = 'on';
                    opts.StartPoint = [1 0.5]; % [Ymin k]
                    opts.Upper = [2 0.1]; % [Ymin k]

                    % Fit model to data.
                    dosesincluding0 = [0;doses{c}];

                    nr = 6;

                    if j == 1
                        [fitresult1{d,c}, gof1(d,c)] = fit( xData, yData, ft, opts );

                        h=subplot(2,nr,d);
                        plot( fitresult1{d,c}, xData, yData );
                        vectcoeff=coeffvalues(fitresult1{d,c});
                        parameters=strcat(['Curve Fit',newline,'Ymax = ',num2str(yMax,'%6.1f'),newline,'\color{red}k = ',num2str(vectcoeff(2),'%6.3f'),newline,'\color{black}R^2 = ',num2str(gof1(d,c).rsquare,'%6.3f')]);
                        ylabel(h,[]);
                        xlabel(h,[]);

                        % assign k values to allkvalues variable which is saved to workspace
                        all_k1{d,c} = vectcoeff(2);
                        allymax1{d,c} = yMax;
                        norm_k = all_k1{d,c}./all_k1{1,c};
                        norm_ymax = allymax1{d,c}./allymax1{1,c};
                        relativegrowthrates1{d,c} = norm_k*norm_ymax;
                        Rsq1{d,c} = (gof1(d,c).rsquare);
                    else
                        [fitresult2{d,c}, gof2(d,c)] = fit( xData, yData, ft, opts );

                        h=subplot(2,nr,d);
                        plot( fitresult2{d,c}, xData, yData );
                        vectcoeff=coeffvalues(fitresult2{d,c});
                        parameters=strcat(['Curve Fit',newline,'Ymax = ',num2str(yMax,'%6.1f'),newline,'\color{red}k = ',num2str(vectcoeff(2),'%6.3f'),newline,'\color{black}R^2 = ',num2str(gof2(d,c).rsquare,'%6.3f')]);
                        ylabel(h,[]);
                        xlabel(h,[]);

                        % assign k values to allkvalues variable which is saved to workspace
                        all_k2{d,c} = vectcoeff(2);
                        allymax2{d,c} = yMax;
                        norm_k = all_k2{d,c}./all_k2{1,c};
                        norm_ymax = allymax2{d,c}./allymax2{1,c};
                        relativegrowthrates2{d,c} = norm_k*norm_ymax;
                        Rsq2{d,c} = (gof2(d,c).rsquare);
                    end

                    %create legend
                    well = string(dosesincluding0(d));
                    firstlegendentry = strcat(['\bf',well{1},'\bf µM']);
                    legend( h, firstlegendentry, parameters, 'Location', 'NorthWest', 'Interpreter', 'tex','color','w');

                    % Set the remaining axes and box properties
                    ax = gca;
                    grid on;
                    ylim([yRange_min yRange_max]);
                    xticks(0:24:144);
                    set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
                        'XGrid','off','Box','on','Color','none');

                end %loop doses
% 
                titletext = append(cellline{a},' ',drug{c},' ',Plate,newline,'f(x) = a*exp(x*k)');

                % Give common xlabel, ylabel and title to your figure
                han=axes(fig,'visible','off');
                han.Title.Visible='on'; han.XLabel.Visible='on'; han.YLabel.Visible='on';
                xlabel(han,'Time [h]','FontWeight','bold','FontSize',11);
                ylabel(han,yaxisname,'FontWeight','bold','FontSize',11);
                title(han,titletext,'FontSize',11);

                filetext = append(date,'_DR_',cellline{a},'_',Plate,'_',drug{c},'_',metric{b});
                saveas(fig, [ filetext, '.svg']);

            end %loop wells

        end %Loop drug

        k1 = cell2table(all_k1,"VariableNames",drug);
        k2 = cell2table(all_k2,"VariableNames",drug);
        Rsq1 = cell2table(Rsq1,"VariableNames",drug);
        Rsq2 = cell2table(Rsq2,"VariableNames",drug);
        RGF1 = cell2table(relativegrowthrates1,"VariableNames",drug);
        RGF2 = cell2table(relativegrowthrates2,"VariableNames",drug);

        writetable(k1,outputfile,'sheet','k1');
        writetable(k2,outputfile,'sheet','k2');
        writetable(Rsq1,outputfile,'sheet','Rsq1');
        writetable(Rsq2,outputfile,'sheet','Rsq2');
        writetable(RGF1,outputfile,'sheet','RGF1');
        writetable(RGF2,outputfile,'sheet','RGF2');

    end %loop metric

end %loop cell lines

end %function