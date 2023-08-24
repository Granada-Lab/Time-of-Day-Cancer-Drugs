function clustergrams_drug_sensitivity(color)

cellline = {'MDAMB468';'HCC1806';'MCF10'};
values = {'GEC50';'GR50';'AOC';'Hill';'GRinf'};
drug  = {'5FU';'Olaparib';'Torin2';'Doxorubicin';'Paclitaxel';'Alisertib';'Cisplatin'};
path = '/Users/carolinector/Nextcloud/Manuscripts/ToDMethods/Data/Fig3/fig3_data/';
file = 'clustergrams_data.xlsx';
sheet = 'ExpFit_CellNr';

rows = [1,2,3]; %rows = cell lines
cols = [1,8,15,22,30]; %cols = sensitivity value

%load data
pathtodata = append(path,file,sheet);
[num] = xlsread(pathtodata);

%% clustergram
v = 1; %sensitivity value

val = values{v};
col = cols(:,v);

for c = 1:3 %collect data for all three cell lines
    row = rows(:,c);
    allvalues(c,:) = num(row,(col:col+6));
end

Y = transpose(allvalues);
cgo = clustergram(Y,'Colormap',color); %@knnimpute,'Standardize','Row'
set(cgo,'RowLabels',drug,'ColumnLabels',cellline,'AnnotColor','k','AnnotPrecision',3);
set(cgo,'ColumnLabelsRotate',90,'RowLabelsRotate',0,'Displayrange',0.3);
addTitle(cgo,val);

%% heatmap

% fig = figure;
% fig.Position = [394,502,389,295];
% fig.OuterPosition = [394,502,389,374];
% fig.InnerPosition = [394,502,389,295];

%% version 1: all celllines, all drugs, single value
% Roundedvalues = round(allvalues,3);
% Y = transpose(Roundedvalues);
% h = heatmap(cellline,drug,Y);
% h.Title = val;
% h.ColorScaling = 'scaled';
% h.FontSize = 14;
% h.CellLabelFormat = '%0.2g';
% s = struct(h);
% s.XAxis.TickLabelRotation = 90;  %vertical orientation
% colormap(color);

clear allvalues

end
