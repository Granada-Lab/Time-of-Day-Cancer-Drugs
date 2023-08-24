% Normalize Excel Data and plot ToD Response Curves

function Plot_Growth_Plots_GRinf

%% Define variables
%create folder+subfolders

metric = 'CellNr';
font = {'Helvetica Neue'};
%drug = {'DMSO';'5FU';'Olaparib';'Torin2';'Doxorubicin';'Paclitaxel';'Alisertib'};
drug = {'5FU';'Olaparib';'Torin2';'Doxorubicin';'Paclitaxel';'Alisertib';'Cisplatin'};
% linestyle = {'--';'-';'-';'-';'-';'-';'-';'-';'-';'-';'-';'-'};
linestyle = {'-';'-';'-';'-';'-';'-';'-';'-';'-';'-';'-';'-'};

%% Load and normalize data

yaxisname = 'Cell Number';
ylabeltext1 = append(yaxisname);

% Load excel file where data is stored
path = 'MDAMB468_doses_at_GRinf_and_GEC50.xlsx';
sheet1 = 'CellNr_1';
sheet2 = 'CellNr_2';
sheet3 = 'CellNr_Cis';

[num1] = xlsread(path,sheet1);
[num2] = xlsread(path,sheet2);
[num3] = xlsread(path,sheet3);

meanx=num1(:,1);

color = lines(numel(drug));
% black = [0,0,0]; 
% color  = [black;colorarray];

for c = 1:7 %loop c drug

    if c < 7
        ysmooth1=smooth(num1(:,c+2));
        ysmooth2=smooth(num2(:,c+2));

        ydata_smooth{c} = mean([ysmooth1,ysmooth2], 2);
        stdev_smooth{c} = std([ysmooth1,ysmooth2],[],2);

        % 1. normalize ydata to timepoint 0
        y_normx0_1_all{c} = ysmooth1./ysmooth1(1,:);
        y_normx0_2_all{c} = ysmooth2./ysmooth2(1,:);
        ydata_normx0{c} = mean([y_normx0_1_all{c},y_normx0_2_all{c}], 2);
        stdev_normx0{c} = std([y_normx0_1_all{c},y_normx0_2_all{c}],[],2);

    else %Cisplatin
        
        ysmooth = smooth(num3(:,3));
        stdsmooth = smooth(num3(:,4));

        ydata_smooth{c} = ysmooth;
        stdev_smooth{c} = stdsmooth;

        % 1. normalize ydata to timepoint 0
        ydata_normx0{c} = ysmooth./ysmooth(1,:);
        stdev_normx0{c} = stdsmooth./ysmooth(1,:);
        xcis = num3(:,1);

    end

end


%% Figure1
fig1 = figure(1);

for j = 1:7

    yraw = cell2mat(ydata_smooth(:,j));
    errorraw = cell2mat(stdev_smooth(:,j));

    if j < 7 
        time = transpose(meanx);
    else 
        time = transpose(xcis);
    end

    y1 = transpose(yraw);
    error1 = transpose(errorraw);

    patch([time fliplr(time)], [(y1-error1)  (fliplr(y1+error1))], color(j,:), 'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold all
    plot(time,y1,'LineWidth',2.5,'LineStyle',linestyle(j,:),'Color',color(j,:)); hold on

    varstoclear1 = {'yraw';'errorraw';'time';'y1';'error1'};
    clear(varstoclear1{:});

end

% Create labels and title

xlabel('Time [h]','FontSize',22,'FontName',string(font));
ylabel(ylabeltext1,'FontSize',22,'FontName',string(font));

title('MDAMB468 Growth at c(GR_{inf}))','FontSize',22,'FontName',string(font),'FontWeight','normal');

% Set the remaining axes and box properties
ax = gca;
grid on;
xticks(0:24:144);
ax.YAxis.Exponent = 3;


set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
    'XGrid','off','Box','on','Color','none','FontSize',22,'FontName',string(font));

% add legend
lgd = legend(drug);
set(lgd,'FontSize',15,'Orientation','vertical','Location','northeast','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');
title(lgd,'Drug','FontWeight','normal','FontName',string(font),'FontSize',16);

%save figure1
saveas(fig1, 'GRinf_noCTRL_MDAMB468_Growth_Curves_Raw_Data.svg');

%% Figure 2
fig2=figure(2);

for k = 1:7

    if k < 7 
        time = transpose(meanx);
    else 
        time = transpose(xcis);
    end

    ynormx0 = cell2mat(ydata_normx0(:,k));
    errorx0 = cell2mat(stdev_normx0(:,k));

    y2 = transpose(ynormx0);
    error2 = transpose(errorx0);

    patch([time fliplr(time)], [(y2-error2)  (fliplr(y2+error2))], color(k,:), 'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold all
    plot(time,y2,'LineWidth',2.5,'LineStyle',linestyle(k,:),'Color',color(k,:)); hold on

    varstoclear2 = {'ynormx0';'errorx0';'y2';'error2'};
    clear(varstoclear2{:});

end

% Create labels and title
ylabeltext2 = append('Normalized ',yaxisname);
xlabel('Time [h]','FontSize',22,'FontName',string(font));
ylabel(ylabeltext2,'FontSize',22,'FontName',string(font));

title('MDAMB468 Growth at c(GR_{inf})','FontSize',22,'FontName',string(font),'FontWeight','normal');

ax = gca;
grid on;
xticks(0:24:144);
%ylim([0 Inf]);
set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
    'XGrid','off','Box','on','Color','none','FontSize',22,'FontName',string(font));

% add legend
lgd = legend(drug);
set(lgd,'FontSize',15,'Orientation','vertical','Location','northeast','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');
title(lgd,'Drug','FontWeight','normal','FontName',string(font),'FontSize',16);

%save figure 2
saveas(fig2, 'GRinf_noCTRL_MDAMB468_Growth_Curves_Normalized_Data.svg');

end