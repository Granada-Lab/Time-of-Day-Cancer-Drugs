% Normalize Excel Data and plot ToD Response Curves 

function Plot_DR_Growth_Plots(folder,date,cellline,drug,doses,v,colorarray)

%% Define variables 
%create folder+subfolders

experiment = str2num(date);
metric = {'CellNr';'Conf'};
font = {'Helvetica Neue'};
linestyle = {'--';'-';'-';'-';'-';'-';'-';'-';'-';'-';'-';'-'};

aa = numel(cellline);
bb = numel(metric);
cc = numel(drug);
dd = (numel(doses{1})+1); %+1 = control

ydata_smooth = cell( 1, dd );
stdev_smooth = cell ( 1, dd );

ydata_normx0 = cell( 1, dd ); 
stdev_normx0 = cell( 1, dd ); 

finalcellcount_1 = cell( dd, cc );
finalcellcount_2 = cell( dd, cc );

ydata_normctrl = cell( 1, dd );
stdev_normctrl = cell( 1, dd );

finalresponse = cell( dd, cc );
sd_finalresponse = cell( dd, cc );

% colorarray = copper(dd); black = [0,0,0]; 
% colorarray(1,:) = []; colorarray = flip(colorarray);
black = [0,0,0]; 
color  = [black;colorarray];

%% Load and normalize data
for b = 1:1 %loop a metric{a}

    close all

    if b==1
        yaxisname = 'Cell Number';
        ylabeltext1 = append(yaxisname);
    else 
        yaxisname = 'Confluency';
        ylabeltext1 = append(yaxisname,' [%]');
    end 

    for a = 2:2 %loop b cell lines
      
        % Load excel file where data is stored
        path = append(folder,date,'_DR_',cellline{a},'.xlsx');
        sheet1 = append(metric{b},'_1');
        if experiment == 2021
            sheet2 = append(metric{b},'_1');
        else 
            sheet2 = append(metric{b},'_2');
        end
    
        [num1] = xlsread(path,sheet1);
        [num2] = xlsread(path,sheet2);
        
        if experiment == 20220317 && b == 1 && a == 3
            num1(1:8,:) = [];
            num2(1:8,:) = [];
        end

        for c = 2:2 %loop c drug
    
            %order: CTRL - lowestdose, dose2, ... highestdose)
            if experiment == 20220204
                yvalues = [v(c),v(c)+5,v(c)+4,v(c)+3,v(c)+2,v(c)+1];
            elseif experiment == 20220317
                yvalues = flip([v(c):7:45,3]);
            elseif experiment == 20221026
                %CONTROL
                if c == 9
                    ctrl_1 = 61; ctrl_2 = 61; %plate1; plate2
                elseif c == 7 || c == 8
                    ctrl_1 = 46; ctrl_2 = 25;
                elseif c == 4 || c == 5 || c == 6
                    ctrl_1 = 3; ctrl_2 = 25;
                else
                    ctrl_1 = 3; ctrl_2 = 3;
                end
                %TREATED
                w = [4,11,18,26,33,40,47,54,62]; %plate2
                yvalues_1 = [ctrl_1,v(c),v(c)+1,v(c)+2,v(c)+3,v(c)+4,v(c)+5,v(c)+6];
                yvalues_2 = [ctrl_2,w(c),w(c)+1,w(c)+2,w(c)+3,w(c)+4,w(c)+5,w(c)+6];
            elseif experiment == 2021
                if c == 1
                    yvalues = [2,11,9,7,5,3,12,10,8,6,4];
                else 
                    yvalues = [2,13,11,9,7,5,3,12,10,8,6];
                end 
            end
    
            %load x-values
            if experiment ~= 2021
                t = 2;
            else 
                t = 1;
            end
            xdata1=num1(:,t); 
            xdata2=num2(:,t); 
            xdata=[xdata1,xdata2];
            meanx = mean(xdata, 2);
    
            for d = 1:dd %loop i concentrations
    
                % load and smooth y-data
                if experiment == 20221026
                    ysmooth1=smooth(num1(:,yvalues_1(d)));
                    ysmooth2=smooth(num2(:,yvalues_2(d)));
                else 
                    ysmooth1=smooth(num1(:,yvalues(d)));
                    ysmooth2=smooth(num2(:,yvalues(d)));
                end

                ydata_smooth{d} = mean([ysmooth1,ysmooth2], 2);
                stdev_smooth{d} = std([ysmooth1,ysmooth2],[],2);

                % 1. normalize ydata to timepoint 0
                y_normx0_1_all{d} = ysmooth1./ysmooth1(1,:);
                y_normx0_2_all{d} = ysmooth2./ysmooth2(1,:);
                ydata_normx0{d} = mean([y_normx0_1_all{d},y_normx0_2_all{d}], 2);
                stdev_normx0{d} = std([y_normx0_1_all{d},y_normx0_2_all{d}],[],2);

                finalcellcount_1{d,c} = median(y_normx0_1_all{d}(end-3:end));
                finalcellcount_2{d,c} = median(y_normx0_2_all{d}(end-3:end));
                
                % 2. normalize ydata_normx0 to control over time
                yval_normctrl_1 = y_normx0_1_all{d}./y_normx0_1_all{1};
                yval_normctrl_2 = y_normx0_2_all{d}./y_normx0_2_all{1};
                ydata_normctrl{d} = mean([yval_normctrl_1,yval_normctrl_2], 2);
                stdev_normctrl{d} = std([yval_normctrl_1,yval_normctrl_2],[],2);
    
                %3. extract final response
                finalresponse{d,c} = median(ydata_normctrl{d}(end-3:end));
                sd_finalresponse{d,c} = median(stdev_normctrl{d}(end-3:end));
    
            end %loop concentration
           
            %% Figure1
            fig1 = figure(1);
    
            y_raw_mat = cell2mat(ydata_smooth);
            sd_raw_mat = cell2mat(stdev_smooth);
    
            for j = 1:dd
    
                yraw = y_raw_mat(:,j);
                errorraw = sd_raw_mat(:,j);
            
                time = transpose(meanx);
                y1 = transpose(yraw);
                error1 = transpose(errorraw);
    
                patch([time fliplr(time)], [(y1-error1)  (fliplr(y1+error1))], color(j,:), 'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold all
                plot(time,y1,'LineWidth',2,'LineStyle',linestyle(j,:),'Color',color(j,:)); hold on
            
            end 
    
            % Create labels and title

            xlabel('Time [h]','FontSize',22,'FontName',string(font));
            ylabel(ylabeltext1,'FontSize',22,'FontName',string(font));
    
            titletext = append(cellline{a});
            title(titletext,'FontSize',22,'FontName',string(font),'FontWeight','normal');
    
            % Set the remaining axes and box properties
            ax = gca;
            grid on;
            xticks(0:24:144);

            if b==2
                yticks(0:20:100); %ylim([0 110]);
            else 
                ax.YAxis.Exponent = 3;
            end

            set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
                'XGrid','off','Box','on','Color','none','FontSize',22,'FontName',string(font));
            
            % add legend
            dosesincluding0 = [0;doses{c}];
            legendentries = num2str(dosesincluding0);
            lgd = legend(legendentries);
            set(lgd,'FontSize',15,'Orientation','vertical','Location','northwest','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');
            
            legendtitle=append(drug{c},' [µM]');
            title(lgd,legendtitle,'FontWeight','normal','FontName',string(font),'FontSize',16);
    
            %save figure1
            if experiment == 2021
                destination = append(date,'_DRs_Cisplatin/DR_2021_',cellline{a},'_Figures_improvedfit/',metric{b},'/');
            else 
                destination = append(folder,'DR_',date,'_',cellline{a},'_Figures_improvedfit/',metric{b},'/');
            end

            filetext1 = append(date,'_DR_',cellline{a},'_',drug{c},'_',metric{b},'_Fig1_RawData');
%     
%             saveas(fig1, [destination filetext1, '.png']);
%             saveas(fig1, [destination filetext1, '.svg']);
    
            %% Figure 2
            fig2=figure(2);
    
            ynormx0_mat = cell2mat(ydata_normx0);
            stdev_normx0_mat = cell2mat(stdev_normx0);
            
            for k = 1:dd

                ynormx0 = ynormx0_mat(:,k);
                errorx0 = stdev_normx0_mat(:,k);
            
                y2 = transpose(ynormx0);
                error2 = transpose(errorx0);

                patch([time fliplr(time)], [(y2-error2)  (fliplr(y2+error2))], color(k,:), 'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold all
                plot(time,y2,'LineWidth',2,'LineStyle',linestyle(k,:),'Color',color(k,:)); hold on

            end 
                
            % Create labels and title
            ylabeltext2 = append('Normalized ',yaxisname);
            xlabel('Time [h]','FontSize',22,'FontName',string(font));
            ylabel(ylabeltext2,'FontSize',22,'FontName',string(font));
    
            title(titletext,'FontSize',22,'FontName',string(font),'FontWeight','normal');
    
            ax = gca;
            grid on;
            xticks(0:24:144);
            %ylim([0 Inf]);
            set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
                'XGrid','off','Box','on','Color','none','FontSize',22,'FontName',string(font));
    
            % add legend
            lgd = legend(legendentries);
            set(lgd,'FontSize',15,'Orientation','vertical','Location','northwest','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');
            
            title(lgd,legendtitle,'FontWeight','normal','FontName',string(font),'FontSize',16);
    
            %save figure 2
            filetext2 = append(date,'_DR_',cellline{a},'_',drug{c},'_',metric{b},'_Fig2_Normx0');
    
%             saveas(fig2, [destination filetext2, '.png']);
%             saveas(fig2, [destination filetext2, '.svg']);
           
    
            %% Figure 3
            fig3=figure(3);

            ydata_normctrl_mat = cell2mat(ydata_normctrl);
            stdev_normctrl_mat = cell2mat(stdev_normctrl);
          
            for l = 2:dd

                ynormctrl = ydata_normctrl_mat(:,l);
                errornormctrl = stdev_normctrl_mat(:,l);
            
                y3 = transpose(ynormctrl);
                error3 = transpose(errornormctrl);

                patch([time fliplr(time)], [(y3-error3)  (fliplr(y3+error3))], color(l,:), 'FaceAlpha',0.1, 'EdgeAlpha',0.1, 'HandleVisibility','off'); hold all
                plot(time,y3,'LineWidth',2,'LineStyle',linestyle(k,:),'Color', color(l,:)); hold on

            end 
                
            % Create labels and title
            ylabeltext3 = append('Relative ',yaxisname);
            xlabel('Time [h]','FontSize',22,'FontName',string(font));
            ylabel(ylabeltext3,'FontSize',22,'FontName',string(font));
    
            title(titletext,'FontSize',22,'FontName',string(font),'FontWeight','normal');
            
            % insert line at y=1 (=control)
            yline(1, 'color', color(1,:), 'LineWidth',2,'LineStyle','--');
           
            ax = gca;
            grid on;
            xticks(0:24:144);
            yticks(0:0.25:1); %ylim([0 1.2]);
            set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
                'XGrid','off','Box','on','Color','none','FontSize',22,'FontName',string(font));
        
            % add legend
            legendentries3 = num2str(doses{c});
            lgd = legend(legendentries3);
            set(lgd,'FontSize',15,'Orientation','vertical','Location','southwest','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');
            
            title(lgd,legendtitle,'FontWeight','normal','FontName',string(font),'FontSize',16);
    
            %save figure 3
            filetext3 = append(date,'_DR_',cellline{a},'_',drug{c},'_',metric{b},'_Fig3_RelCTRL');
            filename3 = strcat(filetext3);
    
%             saveas(fig3, [destination filename3, '.png']);
%             saveas(fig3, [destination filename3, '.svg']);

%             close all

        end %Loop drug
    end %loop cell lines
end %loop metric
end %function