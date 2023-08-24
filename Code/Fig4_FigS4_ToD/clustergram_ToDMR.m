function clustergram_ToDMR(file,celllines_clustergram,channel,drugs,color)

%Carolin Ector, 23.08.2023

%Time-of-Day-Cancer-Drugs Manuscript Fig. 4f

%create clustergram of ToDMR values for all drug-cellline combinations

%input: stored in ToD_analysis_workspace.mat
% file: excel file where ToDMR values are stored
% celllines_clustergram: names of the cell lines included in the clustergram
% channel: either brightfield (confluency) or fluorescent (cell number) channel from live-imaging, in the manuscript we focus on the fluorescent channel 
% drugs: names of drugs the cell lines have been treated with
% color: color for clustergram map

drugs(end,:) = [];

for a = 1:2 %loop a channel

    sheet = append('ToDMR_',channel{a});

    [data] = readmatrix(file,'Sheet',sheet);
    data(9,:) = []; %exclude HCC1937_2
    data(:,[1,end]) = []; %exclude first row (all NaN) and Olaparib (only 2 celllines)

    %% clustergram of ToDMR values

    Y=data';
    cgo = clustergram(Y,'Colormap',color,'ImputeFun', @knnimpute);
    set(cgo,'RowLabels',drugs,'ColumnLabels',celllines_clustergram,'AnnotColor','k','AnnotPrecision',1);
    set(cgo,'ColumnLabelsRotate',90,'RowLabelsRotate',0,'Displayrange',0.45);
    
end %loop a channel

end %function
