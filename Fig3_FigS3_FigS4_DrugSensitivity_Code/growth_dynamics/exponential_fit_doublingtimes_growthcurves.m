function exponential_fit_doublingtimes_growthcurves(file,cellline)

%Carolin Ector, 25.08.2023

%Function ...:
% 1. normalizes growth curves of different cell lines
% 2. fits an exponential function to growth curves and extracts the growth rates
% 3. calculates doubling times
% 4. saves calculated values in an separate excel sheet

%Time-of-Day-Cancer-Drugs Manuscript: Fig 3e + calculation of data for  Fig. 3f

%input: stored in "growth_analysis_workspace.mat"
% file: excel file where growth data is stored
% cellline: names of the cell lines with growth data

%% Define variables
cc = numel(cellline);
channel = {'CellNr';'Conf'};
yaxisnames = {'Cell Number';'Confluency (%)'};
legendtext = {'Cell Number_{t0} = ';'Confluence_{t0} = '};

for a = 1:2 %loop a channel

    %create excel file where data will be stored
    outputfile = append('Growth_analysis_results.xlsx');
    yaxisname = append('Normalized ',yaxisnames{a});

    %load data
    sheet1 = append('Data_',channel{a}); %raw growth data
    sheet2 = append('Error_',channel{a}); %std error (across multiple images per cell line)
    [data] = xlsread(file,sheet1);
    [error] = xlsread(file,sheet2);

    fig1 = figure;
    fig1.Position = [1,1,2560,1361];

    %create empty arrays to store results of the exponential fit
    fitresult = cell(cc,1);
    gof = struct( 'sse', cell(cc,1), 'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

    for c = 1:cc %loop c cell lines

        %account for different imaging intervals
        if ismember(c, [1, 6, 9])
            time = (0:1.5:96)';
        else
            time = (0:2:96)';
        end

        %take data from the first 4 days to make different experiments with different recording lengths comparable
        endrow=size(time,1);

        %load, smooth and normalize data
        rawdata = data(1:endrow,c);
        smootheddata = smooth(rawdata,0.35,'rloess');
        normalizeddata = smootheddata./smootheddata(1,:);

        stdev = error(1:endrow,c);
        smoothedstdev = smooth(stdev,0.35,'rloess');
        normalizedstdev = smoothedstdev./smootheddata(1,:);

        normalizeddata_all{a,c} = normalizeddata;
        stdnormalizeddata_all{a,c} = normalizedstdev;
        time_all{c} = time;

        smootheddata_all{a,c} = smootheddata;
        stdsmootheddata_all{a,c} = smoothedstdev;

        %define initial and final cell number or confluency
        ystart = median(normalizeddata(1:3,:)); %initial
        SE_ystart = median(normalizedstdev(1:3,:)); %corresponding error
        yend = median(normalizeddata(end-3:end,:)); %final
        SE_yend = median(normalizedstdev(end-3:end,:)); %corresponding error

        % calculate doubling times
        growthmetrics{6,c} = 96*log(2)/log(yend/ystart); %doubling time
        growthmetrics{7,c} = (96 * log(2) / (log(yend / ystart))^2) * sqrt((SE_yend / yend)^2 + (SE_ystart / ystart)^2); %corresponding error

        % save initial cell count
        growthmetrics{5,c} = median(smootheddata(1:3,:));

        %% Fit expontial function to growth curve

        [xData, yData] = prepareCurveData(time, normalizeddata);
        equation = 'a*exp(x*k)'; %a=starting point, k=growth rate

        % Set up fittype and options.
        ft = fittype( equation, 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Lower = [0.9 -0.1]; % [a k]
        opts.Robust = 'on';
        opts.StartPoint = [1 0.5]; % [a k]
        opts.Upper = [1.1 0.1]; % [a k]

        % Fit model to data.
        [fitresult{c}, gof(c)] = fit( xData, yData, ft, opts );
        vectcoeff=coeffvalues(fitresult{c});

        % save growth rates with corresponding confidence bounds and Rsq value of the fit
        growthmetrics{4,c} = (gof(c,1).rsquare);
        growthmetrics{1,c} = vectcoeff(2);
        confBounds = confint(fitresult{c});
        growthmetrics{2,c} = confBounds(1,2); %lower confidence bounds growthrate
        growthmetrics{3,c} = confBounds(2,2); %upper confidence bounds growthrate

        if mod(cc, 2) == 0
            nr = cc/2;
        else
            nr = (ceil(cc/2)*2)/2;
        end

        % plot exponential fit
        h=subplot(2,nr,c);
        hold all

        x = transpose(xData);
        y = transpose(yData);
        e = transpose(normalizedstdev);

        patch([x, fliplr(x)], [y - e, fliplr(y + e)], 'b','FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'HandleVisibility', 'off');
        plot( fitresult{c}, xData, yData );

        % add parameters to legend of each subplot
        parameters=strcat([legendtext{a},' = ',num2str(growthmetrics{5,c},'%6.1f'),newline,'\color{black}Doubling Time = ',num2str(growthmetrics{6,c},'%6.1f'),' h',newline,'\color{red}k = ',num2str(growthmetrics{1,c},'%6.3f'),newline,'\color{black}R^2 = ',num2str(growthmetrics{4,c},'%6.3f')]);

        ylabel(h,[]);
        xlabel(h,[]);

        %create legend
        firstlegendentry = cellline{c};
        legend( h, firstlegendentry, parameters, 'Location', 'northoutside', 'Interpreter', 'tex','color','w');

        % Set the remaining axes and box properties
        ax = gca;
        grid on;
        xticks(0:24:120);
        set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
            'XGrid','off','Box','on','Color','none','FontSize',14);

        varstoclear = {'xData';'yData';'rawdata';'smootheddata';'normalizeddata';'stdev';'smoothedstdev';'normalizedstdev';'vectcoeff';'time'};
        clear(varstoclear{:});

    end %loop c cell lines

    % Give common xlabel, ylabel and title to your figure
    han=axes(fig1,'visible','off');
    titletext = append('fitfunc ','[ f(x) = ',equation,' ]',newline,' ');
    han.Title.Visible='on'; han.XLabel.Visible='on'; han.YLabel.Visible='on';
    xlabel(han,'Time (h)','FontWeight','bold','FontSize',16);
    ylabel(han,yaxisname,'FontWeight','bold','FontSize',16);
    title(han,titletext,'FontSize',16);

    hold off

    %save figure
%     filetext = append(fitfunc,'_selected_celllines_manuscript_',metric{a},'_v2');
%     saveas(fig1, [ filetext, '.svg']);

    %save values
    varnames = cellline; %column name
    Value = {'Growth rate';'CI growthrate low';'CI growthrate high';'Rsq exp. Fit';'Initial Value';'Doubling time (h)';'Doubling time Std. Error'}; %row name
    outputsheet = append('Metrics_',channel{a});

    T_Value = cell2table(Value);
    T_growthmetrics = cell2table(growthmetrics,"VariableNames",varnames);

    %merge tables
    T_final = [T_Value,T_growthmetrics];

    %save to excel file
    writetable(T_final,outputfile,'sheet',outputsheet);

end % loop a channel
end % function
