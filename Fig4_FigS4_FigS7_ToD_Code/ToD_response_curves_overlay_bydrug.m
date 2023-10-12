function ToD_response_curves_overlay_bydrug(cellline2,color_all,drug,date_all,mcf10_col)

% Carolin Ector, 23.08.2023

% Function reads relative ToD response data from excel sheet (calculated with function 'ToD_growth_plots_parameter_extraction')
    % and overlays ToD profiles for different celllines drug-by-drug

%Time-of-Day-Cancer-Drugs Manuscript Fig.4e, Fig.S4b

%input: stored in "ToD_overlay_workspace.mat"
% date_all: dates of all ToD experiments
% cellline2: names of the cell lines being analysed
% drug: names of drugs used for time-of-day treatments
% color_drugs: color-coding for the different drugs used in the manuscript
% mcf10_col: column order to load mcf10a data in the drug (d) loop

%Define remaining variables
channel = {'Conf';'CellNr'}; %confluency = brightfield channel, cell number = red fluorescent channel from live-imaging
yaxisnames = {'Confluency';'Cell Number'};
timepoints1 = [0;4;8;16;20;24]; %timepoints of treatment for all drugs except cisplatin
timepoints2 = [0;3;6;9;15;18;21;24]; %timepoints of treatment for cisplatin
cc = numel(cellline2);
dd = numel(drug);

%% Load and normalize data

for a = 1:2 %loop a channel

    sheet = append('responses_',channel{a});

    for d = 1:dd %loop d drug

        fig = figure;
        fig.Position = [1073,574,651,420];

        drugtext = drug{d};

        if d == 9
            aa = 10;
        else
            aa = 1;
        end

        for c = aa:cc %loop c cell lines

            if c < 4
                date = date_all{1};
            elseif c > 3 && c < 7
                date = date_all{2};
            elseif c > 6 && c < 10
                date = date_all{3};
            elseif c >= 10
                date = date_all{4};
            end

            if c ~= 10 && c ~= 2 && d == 8 && c ~= 11
                numbersToCheck = [1,6,8];
                if ismember(c, numbersToCheck)
                    date = date_all{6};
                else
                    date = date_all{5};
                end
            end

            if d == 8 && c ~= 10 && c ~= 2 && c ~= 11
                timepoints = timepoints2;
            else
                timepoints = timepoints1;
            end

            % Load excel file where data is stored
            inputfile = append('ToD_experiments_results/',date,'_ToD_Results_4d_',cellline2{c},'.xlsx');
            [data] = readmatrix(inputfile,'Sheet',sheet);

            if c == 10 || c == 11
                col = mcf10_col(:,d);
            elseif c < 10 && d == 8
                col = 2;
            else
                col = d+1;
            end

            if c == 9 && d == 3 || c == 2 && d == 8
                responses(1:6,:) = 1;
            else
                responses = data(:,col);
            end

            %% spline smoothing

            %Spline interpolation provides a smooth curve that passes through the given data points.
            splinefit = fit(timepoints,responses,'smoothingspline','SmoothingParam',0.7);

            % Evaluate spline fit within the range of timepoints
            xRange = linspace(min(timepoints), max(timepoints), 1000); % Adjust the number of points as desired
            yRange = feval(splinefit, xRange);

            %plot results
            h1 = plot(xRange, yRange);

            %fine-tune plots
            linecolor = cell2mat(color_all(c,:));
            set(h1,'LineWidth',2.5,'color',linecolor);
            hold all

            clear xData
            clear yData
            clear responses

        end %loop c celllines

        % add legend
        legendentries = cellline2(aa:cc,:);
        lgd = legend(legendentries);
        lgdOpt = {'FontSize',16,'Location','northeastoutside', ...
            'FontWeight','normal','EdgeColor','none','Color','#f5f5f5'};
        set(lgd,lgdOpt{:});

        % insert line at y=1 (=control)
        yline(1, 'color', [0.4, 0.4, 0.4], 'LineWidth',2.5,'LineStyle','--','HandleVisibility','off');
        hold off

        mainTextOpt = {'FontSize',22,'FontName','Helvetica Neue','FontWeight','normal'};
        axOpt = {'XLimitMethod','padded','YLimitMethod','padded','linewidth',1.5,'YGrid','on','XGrid','off','Box','on','Color','none','FontSize',22};

        xlabel('Time of the Day [h]',mainTextOpt{:});
        ylabeltext = append('Final response relative to ToD-0 ',newline,'[',yaxisnames{a},']');
        ylabel(ylabeltext,mainTextOpt{:});
        title(drugtext,mainTextOpt{:});

        ax = gca;
        grid on;
        xticks(0:4:32); xlim([-1 25]);
        xticklabels({'0','4','8','12','16','20','24','28','32'});
        set(ax,axOpt{:});

        %save figure
        filename = append('ToD_plots/',channel{a},'/overlay_allcelllines_bydrug/',date,'_ToD_',drug{d},'_',channel{a},'_overlay_allcelllines');
        saveas(fig, [ filename '.svg']);

        clear responses
        clear data
        clear legendentries

    end %loop drug
end %loop a channel
end %function
