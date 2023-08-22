function autocorrelation_calculation(file,rep_bmal,rep_per,col_bmal,col_per,celllines)

%Carolin Ector, 22.08.2023, adapted from https://www.mathworks.com/help/econ/autocorr.html

%%Function returns sample autocorrelation function (ACF) of the univariate time series y, including 95.4% ACF confidence bounds
%%For all time series, the lag 0 autocorrelation acf(1) = 1. 

%input
% file: excel sheet where detrended circadian time-series data is stored
% celllines: names of the cell lines being analysed
% rep_bmal / rep_per: number of replicates per reporter cell line
% col_bmal / col_per: columns in the excel sheet where respective time series from the Bmal1-reporter or Per2-reporter is stored 

%load detrended circadian time-series data
[num] = xlsread(file,'detrended');

%define reporter names, colors and remaining variables
reporter = {'Bmal1';'Per2'};
color = {[0.31,0.68,0.09];[0.85,0.09,0.09]};
cc = numel(celllines); %define number of celllines to process in a loop
autocorr_peak = cell(2,cc); %empty cell arrays to store data
autocorr_period = cell(2,cc);

for a = 1:2 %loop j reporter

    if a == 1 %Bmal1
        cc = numel(col_bmal);
        replicates = rep_bmal;
        column = col_bmal;
    elseif a == 2 %Per2
        cc = numel(col_per);
        replicates = rep_per;
        column = col_per;
    end 

    for c = 1:cc %loop c cell line

        rr = replicates(c,:); %define number of replicated to process in a loop

        for r = 1:rr %loop r replicates

            %load data cellline-by-cellline and replicate-by-replicate 
            col = column(c,:)+r-1; 
            data = num(:,col); 

            fig = figure('Visible','off');

            % calculate autocorrelation within the first 5 days of the recordings (=719 lags)
            [acf,lags,bounds] = autocorr(data,NumLags=719);
            h = stem(lags,acf,'filled');
            set(h,'LineStyle','none','color',color{a},'MarkerSize',4); hold on;
            hold on;
            l1 = line(lags,bounds(1)*ones(length(acf),1));
            l2 = line(lags,bounds(2)*ones(length(acf),1));
            set(l1,'color','b','linewidth',1);
            set(l2,'color','b','linewidth',1);

            %identify second peak with period >= 16h (--> exlude first 96 lags)
            acf2 = acf(96:end,:);
            lags2 = lags(96:end,:);
            [pks,locs]=findpeaks(acf2,lags2,'MinPeakdistance',96,'MinPeakHeight',-0.2,'MinPeakProminence',0.03);
            text(locs+.02,pks,num2str((1:numel(pks))')) % to show peaks in plots by arrow

            hold off

            %save values
            autocorr_peak{r,c} = pks(1);
            autocorr_period{r,c} = locs(1)/6; %get period in hours

            %set graph appearance options
            titletext = append(celllines{c},' ',reporter{a},' ',num2str(r));
            xlabel('Lag at 2^{nd} peak (h)','FontSize',20,'FontName','Helvetica Neue');
            ylabel('Autocorrelation at 2^{nd} peak','FontSize',20,'FontName','Helvetica Neue');

            ax = gca;
            grid on;
            ylim([-1 1]);
            xticks(0:144:839);
            xticklabels({'0','24','48','72','96','120'});

            set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
                'XGrid','off','Box','on','Color','none','FontSize',18,'FontName','Helvetica Neue');

            title(titletext,'FontSize',18,'FontName','Helvetica Neue');

            figurename = append('autocorrelation_',celllines{c},'_',reporter{a},'_',num2str(r));
%             saveas(fig,figurename,'svg');

            variabletoclear = {'data';'pks';'locs';'acf';'bounds';'lags'};
            clear(variabletoclear{:})

        end %loop r replicates

    end %loop c cellline

    clear replicates

    sheet1 = append(reporter{a},'_lag');
    sheet2 = append(reporter{a},'_peak');

    if a == 1
        replicate = {'1';'2';'3';'4';'5';'6';'7';'8';'9'};
    else
        replicate = {'1';'2';'3';'4';'5';'6';'7';'8'};
    end

    t_row = cell2table(replicate);
    varnames = celllines(1:cc,:);

    table_lag = cell2table(autocorr_period,'VariableNames',varnames);
    table_peak = cell2table(autocorr_peak,'VariableNames',varnames);
    t_lag = [t_row,table_lag];
    t_peak = [t_row,table_peak];
    writetable(t_lag,'autocorrelation_results.xlsx','Sheet',sheet1);
    writetable(t_peak,'autocorrelation_results.xlsx','Sheet',sheet2);

    clear autocorr_period
    clear autocorr_peak
    clear replicate

end %loop reporter
end %function
