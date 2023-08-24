function cwt_ridgelength_extraction(celllines,path,rep_bmal,rep_per)

%Carolin Ector, 23.08.2023
%extract ridge lengths from continuous wavelet transform spectra

reporter = {'Bmal1';'Per2'};
outputexcel = 'cwt_ridgelengths_threshold_halfmax.xlsx'; 

%input: stored in cwt_ridgelength_workspace.mat
% file: excel sheet where ridge readout data from continuous wavelet transform is stored
% celllines: names of the cell lines being analysed
% rep_bmal / rep_per: number of replicates per reporter cell line

cc = numel(celllines);

%% load and save data
for c = 1:cc %loop c cellline

    if c > 7
        bb = 1; %no Per2 data for the KO cell lines
    else
        bb = 2;
    end

    for b = 1:bb %loop b reporter

        if b == 1 %Bmal1
            if c == 6 %SY5Y
                rr = rep_bmal(:,c);
                rep = [2:1:3,7:1:rr+5]; %no data for SY5Y BMAL1(5), ridge = 0 h %%TH=halfmax
            else
                rr = rep_bmal(:,c); %number of replicates bmal
                rep = 1:1:rr;
            end
        elseif b == 2 %Per2
            if c == 6 %SY5Y
                rr = rep_bmal(:,c);
                rep = [1:1:3,7:1:rr+3]; %no data for SY5Y BMAL1(5), ridge = 0 h %%TH=halfmax
            else
            rr = rep_per(:,c); %number of replicates per
            rep = 1:1:rr;
            end
        end

        for r = 1:rr %loop r replicates

            replicate = num2str(rep(:,r));

            %oad file
            file = append(celllines{c},'_',reporter{b},'_',replicate,'_ridgeRO.csv'); 
            pathtofile = append(path,file);
            [table_pyboatdata] = readtable(pathtofile);
            pyboatdata = table2array(table_pyboatdata);

            time = pyboatdata(:,1);
            data = pyboatdata(:,2);

            %exclude values above 32 hours
            indices = find(abs(data)>32);
            data(indices) = [];
            time(indices) = [];

            %exclude data after sudden jumps
            ipt = findchangepts(data);
            valuebeforeipt = data(ipt-1,:);
            valuebeafteript = data(ipt,:);
            diff = valuebeafteript - valuebeforeipt;
            rowbeforejump = min(find(data == valuebeforeipt));

            if diff > 1
                time = time(1:rowbeforejump,:);
                data = data(1:rowbeforejump,:);
            end

            %process data
            mintime = min(time);

            %add NaN to missing time points
            if mintime ~= 0
                coltoadd1 = numel(0:0.1666666667:mintime);
                emptycols1(1:coltoadd1,:) = NaN;
                data = [emptycols1;data];
                clear coltoadd1
                clear emptycols1
            end

            length2 = size(data,1);

            if length2 > 720
                data = data(1:720,:);
            elseif length2 < 720
                data(end:720,:) = NaN;
            end

            lengthnonans = rmmissing(data);
            
            if c == 6 
                if b == 1 
                    rowSY5Y = [2,3,7,8,9];
                    row = rowSY5Y(:,r);
                elseif b == 2 
                    rowSY5Y = [1,2,3,7,8];
                    row = rowSY5Y(:,r);
                end
            else 
                row = r;
            end

            length_ridge(row,:) = (size(lengthnonans,1)/6);          
           
            vars = {'data','data1','data2','time'};
            clear(vars{:})

        end

       le = size(length_ridge,1);

        if le < 9 && b == 1
            length_ridge(end+1:9,:) = NaN;
        elseif le < 9 && b == 2
            length_ridge(end+1:9,:) = NaN;
        end
        
        if b == 1
            length_ridge_bmal(:,c) = length_ridge;
        else
            length_ridge_per(:,c) = length_ridge;
        end

        if c == 6
            length_ridge_bmal([1,4,5,6,9],c) = 0;
            length_ridge_per([4,5,6],c) = 0;
        end

        clear length_ridge
        clear rep
    end

    if c > 7
        length_ridge_per(1:9,c) = NaN;
    end
    
    ridgelength(1,c) = mean([length_ridge_bmal(:,c),length_ridge_per(:,c)],'all','omitnan');
    ridgelength(2,c) = mean(length_ridge_bmal(:,c),'all','omitnan');
    ridgelength(3,c) = mean(length_ridge_per(:,c),'all','omitnan');
    ridgelength(4,c) = median([length_ridge_bmal(:,c),length_ridge_per(:,c)],'all','omitnan');
    ridgelength(5,c) = median(length_ridge_bmal(:,c),'all','omitnan');
    ridgelength(6,c) = median(length_ridge_per(:,c),'all','omitnan');
    ridgelength(7,c) = std([length_ridge_bmal(:,c),length_ridge_per(:,c)],[],'all','omitnan');
    ridgelength(8,c) = std(length_ridge_bmal(:,c),[],'all','omitnan');
    ridgelength(9,c) = std(length_ridge_per(:,c),[],'all','omitnan');

end

value = {'mixed_mean';'bmal_mean';'per_mean';'mixed_median';'bmal_median';'per_median';'mixed_stdev';'bmal_stdev';'per_stdev'};
rowtable1 = cell2table(value);
varname1 = celllines;
table1 = array2table(ridgelength,'VariableNames',varname1);
table2 = [rowtable1,table1];
writetable(table2,outputexcel,'sheet','Finalvalues_ridgelength');

for n = 1:2

    replicates = {'1';'2';'3';'4';'5';'6';'7';'8';'9'};
    varname2 = celllines;

    if n == 1 %bmal
        datareplicate = length_ridge_bmal;
        sheet2 = 'Bmal1_ridgelength';
    elseif n == 2 %per
        datareplicate = length_ridge_per;
        sheet2 = 'Per2_ridgelength';
    end

    rowtable2 = cell2table(replicates);
    table3 = array2table(datareplicate,'VariableNames',varname2);
    table4 = [rowtable2,table3];
    writetable(table4,outputexcel,'sheet',sheet2);

end