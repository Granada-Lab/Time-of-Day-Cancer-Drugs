% Normalize Excel Data and plot ToD Response Curves

function Plot_DR_GR50_By_Curve_Fit_Hafner2(folder,date,cellline,drug,doses,input_RC1,input_RC2,dosetextexcel,fitfunc)

%% Define variables
% if strcmp( input_RC1,'k1')
%     fitfunc = 'ExpFit';
% elseif  strcmp( input_RC1,'RGF1')
%     fitfunc = 'LogFit';
% end

metric = {'CellNr';'Conf'};
font = {'Helvetica Neue'};

cc = numel(drug);
aa = numel(cellline);
bb = numel(metric);
ii = (numel(doses{1})+1); %+1 = control

GEC50 = cell ( 6,cc ); %row1=plate1, row2=plate2, row3=mean of plates, row4 =stdev across plates
AUC = cell ( 6,cc ); %see GR50 for organization of cell array
AOC = cell ( 6,cc ); %see GR50 for organization of cell array
GRinf = cell ( 6,cc ); %see GR50 for organization of cell array
HillCoeff = cell ( 6,cc ); %see GR50 for organization of cell array
GR50 = cell ( 6,cc ); %see GR50 for organization of cell array

fitresult = cell( 1, cc ); %this is saved to the workspace in an answer array and it contains the equation and the calculated values
gof = struct( 'sse', cell( 1, cc ), 'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] ); %gof is somehow not saved to the workspace

%% Load and normalize data

for a = 5:5 %loop cellline

    for b = 1:1 %loop metric

        if b == 1
            measurement = '[Cell Number]';
        else
            measurement = '[Confluency]';
        end

        %create excel file where data will be stored
%         filename=append(folder,date,'_Results_',fitfunc,'_',cellline{a},'_',metric{b},'.xlsx');
        filename='/Users/carolinector/Nextcloud/GranadaLab/Users/Carolin/Analysis/IncuCyte/TNBC/DoseResponse_Experiments/2021_DRs_Cisplatin/20210721_DR_ExcelSheets/2021_Results_ExpFit_high_seed_MDAMB468_CellNr.xlsx';
%         filename=append(folder,date,'_Results_',fitfunc,'_',cellline{a},'_',metric{b},'_weighted','.xlsx');

%         filename=append(folder,date,'_Results_GR_CellCount_',cellline{a},'_',metric{b},'.xlsx');
%       
%         filename = append(folder,'Results_ExpFit_Olaparib_CellNr_Figure3.xlsx');

        %Load variables
        [num1] = readtable(filename,'sheet',input_RC1);
        [num2] = readtable(filename,'sheet',input_RC2);
        %load weights
%         [R1] = readmatrix(filename,'Sheet','Rsq1');
%         [R2] = readmatrix(filename,'Sheet','Rsq2');

        experiment = str2num(date);

%         figure;

        for c = 1:cc %loop c drug

            doserange = doses{c};

            close all

            k1 = table2array(num1(:,c)); %plate1
            k2 = table2array(num2(:,c)); %plate2

            for j = 1:ii
                GR1(j,c) = 2^[ k1(j,:)/k1(1,:) ] - 1;
                GR2(j,c) = 2^[ k2(j,:)/k2(1,:) ] - 1;
            end
% 
%             GR1(j,c) = filloutliers(GR1(j,c),"linear");
%             GR1(j,c) = filloutliers(GR2(j,c),"linear");

%             GR1(:,c) = table2array(num1(:,c));
%             GR2(:,c) = table2array(num2(:,c)); 


%             if experiment == 20221026 && b == 2 && fitfunc == "LogFit"
%                 if a == 1
%                     if c == 1
%                         GR1(6,c) = NaN;
%                         GR1(:,c) = fillmissing(GR1(:,c), 'linear');
%                     elseif c == 2
%                         GR2(3,c) = NaN;
%                         GR2(:,c) = fillmissing(GR2(:,c), 'linear');
%                     elseif c == 3
%                         GR1([2,3],c) = NaN;
%                         GR2(3,c) = NaN;
%                         GR1(:,c) = fillmissing(GR1(:,c), 'linear');
%                         GR2(:,c) = fillmissing(GR2(:,c), 'linear');
%                     elseif c == 4
%                         GR1(3,c) = NaN;
%                         GR1(:,c) = fillmissing(GR1(:,c), 'linear');
%                     end
%                 elseif a == 2
%                     if c == 3
%                         GR1([5:8],c) = NaN;
%                         GR2([6,7],c) = NaN;
%                     elseif c == 6
%                         GR1(end,c) = NaN;
%                         GR2(end,c) = NaN;
%                     elseif c == 9
%                         GR2(end,c) = NaN;
%                     end
%                     GR1(:,c) = fillmissing(GR1(:,c), 'linear');
%                     GR2(:,c) = fillmissing(GR2(:,c), 'linear');
%                 end
%             end

            mean_GR = mean([GR1(:,c),GR2(:,c)], 2);
            std_GR = std([GR1(:,c),GR2(:,c)],[], 2);

            %GR => col1=plate1,col2=plate2,col3=mean of plates;
            GR = [GR1(:,c),GR2(:,c),mean_GR];
            GR(1,:) = [];
            std_GR(1,:) = [];

%             Rsq1 = R1(2:end,c);
%             Rsq2 = R2(2:end,c);
%             Rsq = mean([Rsq1,Rsq2], 2);

            if experiment == 20221026
                extradescription = '_shortdoserange';
                if c == 4 || c == 9
                    doserange([1,2],:) = [];
                    GR([1,2],:) = [];
                    std_GR([1,2],:) = [];
%                     Rsq([1,2],:) = [];
                elseif c == 8
                    doserange(end,:) = [];
                    GR(end,:) = [];
                    std_GR(end,:) = [];
%                     Rsq(end,:) = [];
                else
                    doserange(1,:) = [];
                    GR(1,:) = [];
                    std_GR(1,:) = [];
%                     Rsq(1,:) = [];
                end
            else
                extradescription = '';
            end

%             Rsq = max(Rsq,0); %replace negative value by 0 
             
            maxdose = max(doserange);
            mindose = min(doserange);

            titletext = append(cellline{a}, ' ',drug{c});

            %controldose = (doses{c}(1,:)/100);
            %Dose = [controldose;doses{c}];
            Dose = [doserange];

            rowinexcel = [1,2,6]; %for k = 1 or 2 or 3

            for k = 1:1 %loop plates & combined

                if k == 1
                    extradescription = '30';
                elseif k == 2
                    extradescription = '60';
                else
                    extradescription = 'mean3060';
                end

                %clear fitresult
                %clear gof

                Response = GR(:,k);
                %
                %                 if experiment == 20221026
                %                     if c == 4 || c == 9
                %                         Response([1,2],:) = [];
                %                     elseif c == 8
                %                         Response(end,:) = [];
                %                     else
                %                         Response(1,:) = [];
                %                     end
                %                 end

                %                 Rsq_raw = Rsq(:,k);
                %                 Weights = max(Rsq_raw, 0); %set negative values to 0

                [xData, yData] = prepareCurveData(Dose, Response);

                % Set up fittype and options.
                ft = fittype( 'GRinf+((1-GRinf)/(1+(c/GEC50)^hGR))', 'independent', 'c', 'dependent', 'GR' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%                 opts.Weights = Rsq;
                opts.Lower = [0.00001 -1  0]; % [GEC50,GRinf,hGR]
                opts.Robust = 'Off';
                opts.StartPoint = [0.5 0.2 1]; % [GEC50,GRinf,hGR]
                opts.Upper = [1000 1 10]; % [GEC50,GRinf,hGR]

                % Fit model to data
                [fitresult{k,c}, gof(k,c)] = fit( xData, yData, ft, opts);
                vectcoeff=coeffvalues(fitresult{k,c});

                row = rowinexcel(:,k);

                test1 = isempty(vectcoeff(1));
                if test1 == 1
                    GEC50{row,c} = 'NaN';
                else
                    GEC50{row,c} = vectcoeff(1);
                end
                AUC{row,c} = trapz(Dose, Response);
                AOC{row,c} = trapz(Dose, (1-Response));
                GRinf{row,c} = vectcoeff(2);
                HillCoeff{row,c} = vectcoeff(3);

                xpoints=logspace(log10(Dose(1,:)),log10(Dose(end,:)),1000);
                ypoints=logspace(log10(Response(1,:)),log10(Response(end,:)),1000);

                %maxresponse = max(Response);
                %minresponse = min(Response);

                if Response(end,:) > 0.502
                    %                     y_around50 = [];
                    %                     closestIndex = [];
                    %                     y_GR50 = [];
                    inf2 = '+';
                    GR50{row,c} = '+infty';
                elseif Response(1,:) < 0.502
                    %                     y_around50 = [];
                    %                     closestIndex = [];
                    %                     y_GR50 = [];
                    inf2 = '-';
                    GR50{row,c} = '-infty';
                else
                    length_responses = size((Response),1);
                    length_uniqueresponses = size(unique(Response),1);
                    diff = length_responses-length_uniqueresponses;
                    if diff ~= 0
                        [C,ia,~] = unique(Response,'rows');
                        for v = 1:length(ia)
                            row_ia = ia(v,:);
                            dose_short(v,:) = Dose(row_ia,:);
                            resp_short(v,:) = Response(row_ia,:);
                        end
                        GR50{row,c} = exp(interp1(resp_short, log(dose_short), 0.5));
                    else
                        GR50{row,c} = exp(interp1(Response, log(Dose), 0.5));
                    end
                    %                     %find y-value closest to y=0.5 (GR50)
                    %                         y_around50 = ypoints(ypoints >= 0.498 & ypoints <= 0.502);
                    %                         [~,closestIndex] = (min(abs(y_around50 - 0.5)));
                    %                         y_GR50 = y_around50(closestIndex);
                    %                         index = find(ypoints==y_GR50);
                    %                         x_GR50 = xpoints(index);
                    %                         %text(x_GR50,y_GR50,[' \leftarrow ' sprintf('GR_{50}=%0.2gµM',x_GR50)],'FontSize',16);
                    %                         GR50{row,c} = x_GR50;
                end
            end

            %evaluate if mean is calculatable by checking for numeric values
            test3_1 = isnumeric(GEC50{1,c}); % test3 = 0 if non-numeric and test3 = 1 for numeric
            test3_2 = isnumeric(GEC50{2,c});
            if test3_1 == 0 || test3_2 == 0 %one or both of the values are non-numeric
                GEC50{3,c} = 'NaN';
                GEC50{4,c} = 'NaN';
            else %both of the values are numeric
                GEC50{3,c} = mean([GEC50{1,c},GEC50{2,c}],2);
                GEC50{4,c} = std([GEC50{1,c},GEC50{2,c}],[],2);
                if GEC50{3,c} > maxdose
                    GEC50{3,c} = '+infty';
                    GEC50{4,c} = 'NaN';
                elseif GEC50{3,c} < mindose
                    GEC50{3,c} = '-infty';
                    GEC50{4,c} = 'NaN';
                end
            end

            test3_3 = isnumeric(GEC50{6,c});
            if test3_3 == 0
                GEC50{6,c} = 'NaN';
            else
                if GEC50{6,c} > maxdose
                    GEC50{6,c} = '+infty';
                elseif GEC50{6,c} < mindose
                    GEC50{6,c} = '-infty';
                end
            end

            AUC{3,c} = mean([AUC{1,c},AUC{2,c}],2);
            AUC{4,c} = std([AUC{1,c},AUC{2,c}],[],2);
            AOC{3,c} = mean([AOC{1,c},AOC{2,c}],2);
            AOC{4,c} = std([AOC{1,c},AOC{2,c}],[],2);
            GRinf{3,c} = mean([GRinf{1,c},GRinf{2,c}],2);
            GRinf{4,c} = std([GRinf{1,c},GRinf{2,c}],[],2);
            HillCoeff{3,c} = mean([HillCoeff{1,c},HillCoeff{2,c}],2);
            HillCoeff{4,c} = std([HillCoeff{1,c},HillCoeff{2,c}],[],2);

            %evaluate if mean is calculatable by checking for numeric values
            test4_1 = isnumeric(GR50{1,c}); % see above for explanation
            test4_2 = isnumeric(GR50{2,c});
            if test4_1 == 1 || test4_2 == 1 %both of the values are numeric
                GR50{3,c} = mean([GR50{1,c},GR50{2,c}],2);
                GR50{4,c} = std([GR50{1,c},GR50{2,c}],[],2);
            end

            error = std_GR;

            fig = figure;%('Visible','off');
            h = plot(fitresult{1,c}); hold all;
            set(h,'lineWidth',1,'Color','r','HandleVisibility','off');
            hold all;
            errorbar(xData,yData,error,'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'LineWidth',1,'Color',[0 0 0],'HandleVisibility','off');

            %plot(xData(:,2),yData(:,2),'o','MarkerSize',7,'MarkerFaceColor',[0 0 0],'LineWidth',1,'Color',[0 0 0],'HandleVisibility','off');

            %dummyline for legend
            line(xData, yData, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none'); %GEC50
            line(xData, yData, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none'); %GR50
            line(xData, yData, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none'); %GRinf
            line(xData, yData, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none'); %Hill
            line(xData, yData, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none'); %AOC
            hold off

            test5_1 = isnumeric(GEC50{3,c});
%             test5_1 = isnumeric(GEC50{6,c});
            %test5_1 = isnumeric(GEC50{2,c});

            if test5_1 == 0
                lgd_GEC50=strcat('GEC_{50} = NaN');
            else
                lgd_GEC50=strcat(['GEC_{50} = ',num2str(GEC50{3,c},'%6.3f'),' ',char(177),' ',num2str(GEC50{4,c},' %6.3f µM')]);
%                 lgd_GEC50=strcat(['GEC_{50} = ',num2str(GEC50{6,c},'%6.3f µM')]);
                %lgd_GEC50=strcat(['GEC_{50} = ',num2str(GEC50{2,c},'%6.3f µM')]);
            end

            test5_2 = isnumeric(GR50{3,c});
%             test5_2 = isnumeric(GR50{6,c});
            %test5_2 = isnumeric(GR50{2,c});
            if test5_2 == 0
                lgd_GR50=strcat('GR_{50} = ',inf2,'\infty');
            else
                lgd_GR50=strcat(['GR_{50} = ',num2str(GR50{3,c},'%6.3f'),' ',char(177),' ',num2str(GR50{4,c},' %6.3f µM')]);
%                 lgd_GR50=strcat(['GR_{50} = ',num2str(GR50{6,c},'%6.3f µM')]);
%                 lgd_GR50=strcat(['GR_{50} = ',num2str(GR50{2,c},'%6.3f µM')]);
            end
            lgd_GRinf=strcat(['GR_{inf} = ',num2str(GRinf{3,c},'%6.3f'),' ',char(177),' ',num2str(GRinf{4,c},' %6.3f µM')]);
            lgd_Hill=strcat(['Hill = ',num2str(HillCoeff{3,c},'%6.3f'),' ',char(177),' ',num2str(HillCoeff{4,c},' %6.3f µM')]);
            lgd_AOC=strcat(['AOC = ',num2str(AOC{3,c},'%6.3f'),' ',char(177),' ',num2str(AOC{4,c},' %6.3f µM')]);

%             lgd_GRinf=strcat(['GR_{inf} = ',num2str(GRinf{6,c},'%6.3f')]);
%             lgd_Hill=strcat(['Hill = ',num2str(HillCoeff{6,c},'%6.3f')]);
%             lgd_AOC=strcat(['AOC = ',num2str(AOC{6,c},'%6.3f')]);
            %             lgd_GRinf=strcat(['GR_{inf} = ',num2str(GRinf{2,c},'%6.3f')]);
            %             lgd_Hill=strcat(['Hill = ',num2str(HillCoeff{2,c},'%6.3f')]);
            %             lgd_AOC=strcat(['AOC = ',num2str(AOC{2,c},'%6.3f')]);

            legend(lgd_GEC50,lgd_GR50,lgd_Hill,lgd_GRinf,lgd_AOC,'Interpreter','tex','color','none','EdgeColor','none','FontSize',14);

            %label axes
            title(titletext,'FontSize',18,'FontName',string(font),'FontWeight','normal');
            xlabel('Dose [µM]','FontSize',18,'FontName','Helvetica Neue');
            ylabeltext = append('GR ',measurement);
            ylabel(ylabeltext,'FontSize',18,'FontName','Helvetica Neue');

            ax = gca;

            set(ax,'XScale','log','XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
                'XGrid','off','Box','on','Color','none','FontSize',18,'FontName',string(font));

            xticks(doserange);
            xticklabels(string(doserange));

            %save figure
%             filetext = append(date,'_DR_',cellline{a},'_',drug{c},'_',metric{b},'_Fig6_GR50_',fitfunc,'_noweights_averagevalues',extradescription);
            filetext = append(date,'_DR_',cellline{a},'_',drug{c},'_',metric{b},'_Fig6_GR50_',fitfunc,'_noweights_averagevalues',extradescription);
            subfolder = append('DR_',date,'_',cellline{a},'_Figures_improvedfit/',metric{b},'/');
%             subfolder = append('DR_',date,'_',cellline{a},'_Figures_GRcellcount/',metric{b},'/');
%             if experiment == 2021
%                 destination = append(date,'_DRs_Cisplatin/DR_2021_',cellline{a},'_Figures_GRcellcount/',metric{b},'/');
%             else
%                 destination = append(folder,subfolder);
%         end
%             destination = append(folder,subfolder);

%             saveas(fig, [ destination filetext, '.png']);
%             saveas(fig, [ destination filetext, '.svg']);
            saveas(fig, [ filetext, '.svg']);
%             saveas(fig, [ destination filetext, '.fig']);

%             clear Rsq1
%             clear Rsq2
%             clear Rsq

        end %loop drug

%         hold off

%         varnames = drug;
% 
%         t_GR1 = array2table(GR1,"VariableNames",varnames,"RowNames",dosetextexcel);
%         t_GR2 = array2table(GR2,"VariableNames",varnames,"RowNames",dosetextexcel);
%         writetable(t_GR1,filename,'sheet','GR1');
%         writetable(t_GR2,filename,'sheet','GR2');
% 
%         %finalvalues = {'Plate1','Plate2','Mean','StdDev',' ','Combined'};
%         t_GEC50 = cell2table(GEC50,"VariableNames",varnames);
%         t_AUC = cell2table(AUC,"VariableNames",varnames);
%         t_AOC = cell2table(AOC,"VariableNames",varnames);
%         t_GRinf = cell2table(GRinf,"VariableNames",varnames);
%         t_HillCoeff = cell2table(HillCoeff,"VariableNames",varnames);
%         t_GR50 = cell2table(GR50,"VariableNames",varnames);
%         writetable(t_GEC50,filename,'sheet','GEC50');
%         writetable(t_AUC,filename,'sheet','AUC');
%         writetable(t_AOC,filename,'sheet','AOC');
%         writetable(t_GRinf,filename,'sheet','GRinf');
%         writetable(t_HillCoeff,filename,'sheet','HillCoeff');
%         writetable(t_GR50,filename,'sheet','GR50');

        clear GEC50;
        clear AUC;
        clear AOC;
        clear GRinf;
        clear HillCoeff;
        clear GR50;

    end %loop metric
end %Loop cellline
end %loop function