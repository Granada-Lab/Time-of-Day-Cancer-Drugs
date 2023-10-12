function ToD_response_curves_overlay_bycellline(date_all,cellline2,drug,color_drugs,mcf10_col)

% Carolin Ector, 23.08.2023

% Function reads relative ToD response data from excel sheet (calculated with function 'ToD_growth_plots_parameter_extraction')
    % and overlays ToD profiles for different drugs cellline-by-cellline

%Time-of-Day-Cancer-Drugs Manuscript Fig.4d, Fig.S4a

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

    for c = 11:cc %loop c cell lines

        if c == 2
            dd = 7;
        elseif c == 10 || c == 11
            dd = 9;
        else
            dd = 8;
        end

        fig = figure;
        fig.Position = [1073,574,639,420];
        titletext = cellline2{c};

        for d = 1:dd %loop d drug

            if c < 4 && d < 8
                date = date_all{1};
            elseif c > 3 && c < 7 && d < 8
                date = date_all{2};
            elseif c > 6 && c < 10 && d < 8
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

            if c == 2
                dd = 7;
            elseif c >= 10
                dd = 9;
            else
                dd = 8;
            end

            if d == 8 && c ~= 10 && c ~= 2 && c ~= 11
                timepoints = timepoints2;
            else
                timepoints = timepoints1;
            end

            % Load excel file where data is stored
            sheet = append('responses_',channel{a});
            inputfile = append('ToD_experiments_results/',date,'_ToD_Results_4d_',cellline2{c},'.xlsx');
            [data] = readmatrix(inputfile,'Sheet',sheet);

            if c == 10 || c == 11
                col = mcf10_col(:,d);
            elseif c < 10 && d == 8
                col = 2;
            else
                col = d+1;
            end

            if c == 9 && d == 3
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
            h1 = plot(xRange, yRange);
            legendentries = drug(1:dd,:);

            %% fine-tune plots
            set(h1,'LineWidth',3,'color',color_drugs(d,:));

            hold all

            clear xData
            clear yData
            clear responses

        end %loop d drug

        % add legend
        lgd = legend(legendentries);
        lgdOpt = {'FontSize',16,'Location','northeastoutside', ...
            'FontWeight','normal','EdgeColor','none','Color','#f5f5f5'};
        set(lgd,lgdOpt{:});

        % insert line at y=1 (=control)
        yline(1, 'color', [0.4, 0.4, 0.4], 'LineWidth',2.5,'LineStyle','--','HandleVisibility','off');
        hold off

        mainTextOpt = {'FontSize',22,'FontName','Helvetica Neue','FontWeight','normal'};
        axOpt = {'XLimitMethod','padded','YLimitMethod','padded','linewidth',1.5,'YGrid','on','XGrid','off','Box','on','Color','none','FontSize',22};

        xlabel('Time of Day (h)',mainTextOpt{:});
        ylabeltext = append('Relative ',yaxisnames{a});
        ylabel(ylabeltext,mainTextOpt{:});
        title(titletext,mainTextOpt{:});

        ax = gca;
        grid on;
        xticks(0:4:32); xlim([-1 25]);
        xticklabels({'0','4','8','12','16','20','24','28','32'});
        set(ax,axOpt{:});

        %% save figure
        filename = append('ToD_plots/',channel{a},'/overlay_alldrugs_bycellline/',date,'_ToD_',cellline2{c},'_',channel{a},'_overlay_alldrugs');
        saveas(fig, [ filename '.svg']);

        clear data
        clear legendentries

    end %loop c cell lines
end %loop a channel
end %function
