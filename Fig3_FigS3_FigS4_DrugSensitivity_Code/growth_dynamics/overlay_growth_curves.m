function overlay_growth_curves(file,cellline,color_all)

%Carolin Ector, 25.08.2023

%function plots growth curves as a function of the recording time:
%1. overlays growth curves from the two different imaging channels
    %(CellNr = fluorescent channel, nuclei counts; Conf = brightfield channel, confluency in %)
%2. overlays normalized growth cruved of different cell lines (normalized to timepoint 0) 

%Time-of-Day-Cancer-Drugs Manuscript Fig.3c,d and Fig.S2a

%input: stored in "growth_analysis_workspace.mat"
% file: excel file where growth data is stored
% cellline: names of the cell lines with growth data
% color_all: color array containing a distinct color for each cell line

channel = {'CellNr';'Conf'};
cc = numel(cellline);

for c = 1:cc %loop c cell lines

    for a = 1:2 %loop a channel

        sheet1 = append('Data_',channel{a});
        sheet2 = append('Error_',channel{a});
        [data] = readmatrix(file,'sheet',sheet1);
        [stderr] = readmatrix(file,'sheet',sheet2);

        %account for different imaging intervals
        if ismember(c, [1, 6, 9])
            time = (0:1.5:96)';
        else
            time = (0:2:96)';
        end

        time_all{:,c} = time;

        %take data from the first 4 days to make different experiments with different recording lengths comparable
        endrow=size(time,1);

        %load, smooth and normalize data
        rawdata = data(1:endrow,c);
        smootheddata = smooth(rawdata,0.35,'rloess');
        smootheddata_all{a,c} = smootheddata;
        normalizeddata_all{a,c} = smootheddata./smootheddata(1,:);

        stderr = stderr(1:endrow,c);
        smoothederror = smooth(stderr,0.35,'rloess');
        smoothedstderr_all{a,c} = smoothederror;
        normalizedstderr_all{a,c} = smoothederror./smootheddata(1,:);

        clear smootheddata
        clear smoothederror
        clear rawdata
        clear time

    end %loop a channel

    fig = figure;%('visible','off');

    %add standard deviation as shaded area to plot
    x = transpose(cell2mat(time_all(:,c)));

    %load data cell number
    y1 = transpose(cell2mat(smootheddata_all(1,c)));
    error1 = transpose(cell2mat(smoothedstderr_all(1,c)));

    %load data confluency
    y2 = transpose(cell2mat(smootheddata_all(2,c)));
    error2 = transpose(cell2mat(smoothedstderr_all(2,c)));

    hold all
    yyaxis left
    patch([x, fliplr(x)], [y1 - error1, fliplr(y1 + error1)],'r','FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'HandleVisibility', 'off');
    plot(x,y1,'Color','r','LineWidth',2);
    ylabel('Cell Number');

    yyaxis right
    plot(x,y2,'Color','b','LineWidth',2);
    patch([x, fliplr(x)], [y2 - error2, fliplr(y2 + error2)],'b','FaceAlpha', 0.1, 'EdgeAlpha', 0.1, 'HandleVisibility', 'off');
    ylabel('Confluency (%)');
    yyaxis right

    xlabel('Time (h)');
    title(cellline{c});

    hold off

    % Set the remaining axes and box properties
    ax = gca;
    grid on;
    xticks(0:24:120);
    set(ax,'XLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
        'XGrid','off','Box','on','Color','w','FontSize',16);

    %save figure
    filetext1 = append('overlay_growthcurves_cellnr_conf_',cellline{c});
    %     saveas(fig, [ filetext1, '.svg']);

    varstoclear = {'x';'y1';'error1';'y2';'error2'};
    clear(varstoclear{:});

end %loop c cell lines

for a2 = 1:2 %loop2 channel

    fig = figure;
    colorder = [10;1;2;3;4;5;6;7;8;9];
    cellline2 = cellline(colorder,:);

    for c2 = 1:cc %loop2 cellline

        %column where data for the cellline is stored in excel sheet
        col = colorder(c2,:);

        x2 = transpose(cell2mat(time_all(1,col)));
        y = transpose(cell2mat(normalizeddata_all(a2,col)));
        error = transpose(cell2mat(normalizedstderr_all(a2,col)));

        color = str2num(cell2mat(color_all(col,:)));

        hold all;
        patch([x2, fliplr(x2)], [y - error, fliplr(y + error)],color,'EdgeColor', 'none','FaceAlpha', 0.1, 'HandleVisibility', 'off');
        plot(x2,y,'Color',color,'LineWidth',3);

    end %loop2 cellline

    ylabeltext = append('Normalized ',channel{a2});
    ylabel(ylabeltext);
    xlabel('Time (h)');

    %create legend
    legend(cellline2, 'Location', 'NorthWest', 'Interpreter', 'tex','color','w');

    % Set the remaining axes and box properties
    ax = gca;
    grid on;
    xticks(0:24:120);
    set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
        'XGrid','off','Box','on','Color','w','FontSize',20);

    %save figure
    filetext2 = append('overlay_normalized_growth_curves_all_celllines');
    % saveas(fig2, [ filetext2, '.svg']);

    varstoclear = {'x2';'y';'error'};
    clear(varstoclear{:});

end %loop2 channel
end %function