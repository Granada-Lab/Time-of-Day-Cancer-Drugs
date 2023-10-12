function DR_growth_plots_GRinf(cellline,channel,date_all,drug,color_drugs)

%Carolin Ector, 27.09.2023

%Function plots growth curves of cell lines treated at their respective closest and smallest dose which elicits a GRinf (maximal) drug response.

%Time-of-Day-Cancer-Drugs Manuscript Fig. 3j

%input: stored in 'workspace_growth_GRinf.mat'
% date_all: dates of the experiments being analyzed
% channel: channel being analyzed from live-imaging: confluency = brightfield channel, cell number = red fluorescent channel
% cellline: names of the cell lines being analysed
% drug: names of drugs used for drug treatments of the dose-response experiments
% color_drugs: colors for the different drugs

cc = numel(cellline);
yaxisnames = {'Normalized Cell Number';'Normalized Confluency'};
font = {'Helvetica Neue'};

for c = 1:cc %loop c cell line

    if c == 3
        d1 = 1;
        dd = numel(drug)-1;
    else
        d1 = 1;
        dd = numel(drug);
    end

    if c == 0
        a1 = 2;
    else
        a1 = 1;
    end

    for a = a1:numel(channel) %loop a channel

        figure;

        for d = d1:dd

            if d == 8
                date = date_all{4};
            else
                if c < 3
                    colorder = [1,3,6,4,5,7,7,1,2];
                    if d == 6 || d == 8
                        date = date_all{1};
                    else
                        date = date_all{2};
                    end
                elseif c == 3
                    colorder = [1,3,4,5,6,7,8,1];
                    date = date_all{1};
                elseif c == 4
                    colorder = [1,3,6,4,5,6,7,1,2];
                    date = date_all{2};
                elseif c == 5 || c == 6
                    colorder = [1,3,9,6,7,4,8,1,2];
                    date = date_all{3};
                end
            end

            %load normalized growth curves from inputfile
            inputfile=append('DR_results/',date,'_Growth_GRinf_',cellline{c},'.xlsx');
            inputsheet1 = append('growth_GRinf_',channel{a});
            inputsheet2 = append('growth_GRinf_std_',channel{a});

            [data] = readmatrix(inputfile,'sheet',inputsheet1);
            [std] = readmatrix(inputfile,'sheet',inputsheet2);

            col = colorder(:,d)+1;

            x = transpose(data(:,1));
            y = transpose(data(:,col));
            err = transpose(std(:,col));

            if c == 4 && col == 7 %no data for Alpelisib and Adavosertib
                y(:,1:end) = 1;
                err(:,1:end) = 0;
            end

            h1 = patch([x fliplr(x)], [(y-err)  (fliplr(y+err))], color_drugs(d,:), 'FaceAlpha',0.1, 'EdgeAlpha',0, 'HandleVisibility','off'); hold all
            plot(x,y,'LineWidth',2,'Color',color_drugs(d,:));

            varstoclear = {'data','std','x','y','err','date'};
            clear(varstoclear{:})

        end %loop d drug

        hold off

        % Create labels and title
        xlabel('Time (h)','FontSize',22,'FontName',string(font));
        ylabel(yaxisnames{a},'FontSize',22,'FontName',string(font));
        title(cellline{c},'FontSize',22,'FontName',string(font),'FontWeight','normal');

        % Set the remaining axes and box properties
        ax = gca;
        grid on;
        xticks(0:24:144);

        set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
            'XGrid','off','Box','on','Color','none','FontSize',22,'FontName',string(font));

        % add legend
        lgd = legend(drug);
        set(lgd,'FontSize',15,'Orientation','vertical','Location','northeast','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');
        title(lgd,'Drug','FontWeight','normal','FontName',string(font),'FontSize',16);

        %save figure
        filetext = append('DR_plots/','_DR_',cellline{c},'_',channel{a},'_Fig7_Growth_GRinf.svg');
        saveas(h1, filetext);

    end %loop a channel
end %loop c cellline
end %function