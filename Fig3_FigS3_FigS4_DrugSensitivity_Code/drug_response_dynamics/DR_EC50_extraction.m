function [parameters] = DR_EC50_extraction(a,c,date,cellline,drug,doses,finalresponse,sd_finalresponse,parameters)

%Carolin Ector, 28.09.2023
%Adapted from Ritchie Smith (2023). doseResponse (https://www.mathworks.com/matlabcentral/fileexchange/33604-doseresponse),MATLAB Central File Exchange.

%Time-of-Day-Cancer-Drugs Manuscript: calculates needed input variables for generation of Fig.3n

%Function computes the EC50 of a dose/response relationship given two vectors describing the doses and responses.
%A semilog graph is plotted illustrating the relationship, upon which the mean and standard error of the response to each dose level is
%plotted along with the fitted Hill Equation sigmoid.Requires nlinfit from statistics toolbox.

%input: partially stored in '[date]_DR_workspace.mat',remaining variables are defined in 'DR_analysis_pipeline.m' script
% a: loop a channel being analyzed from live-imaging: 1 = cell number (red fluorescent channel), 2 = confluency (brightfield channel)
% c: loop c celllines
% date: date of the experiment being analyzed
% channel: channel being analyzed from live-imaging (see above)
% cellline: names of the cell lines being analysed
% drug: names of drugs used for drug treatments of the dose-response experiments
% doses: doses of each drug administered to the cells
% finalresponse: final nucleus counts (relative to control) for each drug dose
% sd_finalresponse: standard deviation of final nucleus counts (relative to control) for each drug dose
% parameters: array to store final drug sensitvity parameters

%output:
% EC50 value stored in parameters(6,:), not: parameters(1:5,:) = drug sensitivty values computed from the Hafner approach

%Define remaining variables
dd = numel(drug);
channel = {'CellNr';'Conf'};
yaxisnames = {'Cell Number';'Confluency'};
experiment = str2num(date);

for d = 1:dd %loop d drug

    Doserange = doses{d};

    if experiment == 2021
        factor = 2;
    else
        factor = 4;
    end

    xaxislabels = [Doserange(1,:)/factor;Doserange;Doserange(end,:)*factor];
    xaxislabels = xaxislabels(1:2:end);
    ylabeltext = append('Relative ',yaxisnames{a});

    %Convert cell array to matrix
    Response = cell2mat(finalresponse(:,d));
    Error = cell2mat(sd_finalresponse(:,d));

    %exclude controls from sigmoidal response curve fit
    Response(1,:) = [];
    Error(1,:) = [];

    if experiment == 20221026 %additional doses tested compared to the other experiments
        if d == 4 || d == 9
            Doserange([1,2],:) = [];
            Response([1,2],:) = [];
            Error([1,2],:) = [];
        elseif d == 8
            Doserange(end,:) = [];
            Response(end,:) = [];
            Error(end,:) = [];
        else
            Doserange(1,:) = [];
            Response(1,:) = [];
            Error(1,:) = [];
        end
    else
    end

    %Excluding ohighest dose and response for Doxorubicin
    %(Doxorubicin-induced red signal enhancement due to drug autofluorescence)
    if strcmp(drug{d}, 'Doxorubicin')
        Response(end,:) = [];
        Doserange(end,:) = [];
        Error(end,:) = [];
    end

    %Define sigmoidal model function
    sigmoid=@(beta,x)beta(1)+(beta(2)-beta(1))./(1+(x/beta(3)).^beta(4));

    %calculate some rough guesses for initial parameters
    minResponse=min(Response); %beta(1)
    maxResponse=max(Response); %beta(2)
    minDose=min(Doserange);
    maxDose=max(Doserange);
    midDose=mean([minDose maxDose]); %beta(3)

    % Define lower and upper bounds for the parameters (constraints)
    lowerBounds = [0, 1, minDose/100, 0.5]; %previously beta(4)=1
    upperBounds = [1, 1, maxDose*100, 10];

    % Define initial parameter guesses
    initialGuess = [minResponse, maxResponse, midDose, 5];

    % Set options for lsqnonlin (you can customize as needed)
    options = optimset('Display', 'off'); % Turn off display for quieter optimization
    [coeffs] = lsqnonlin(@(beta) Response - sigmoid(beta, Doserange), initialGuess, lowerBounds, upperBounds, options);

    % Extract the fitted parameters
    parameters(6,d) = coeffs(3);

    %plot the fitted sigmoid
    fig = figure('Visible','off');

    xpoints=logspace(log10(minDose/factor),log10(maxDose*factor),1000);
    semilogx(xpoints,sigmoid(coeffs,xpoints),'Color',[0 0 0],'LineWidth',1,'HandleVisibility','off')
    hold on

    %dummyline for legend
    line(Doserange, Response, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
    lgd_EC50=strcat(['EC_{50} = ',num2str(coeffs(3),'%6.3f')]);
    legend(lgd_EC50,'Interpreter','tex','color','none','EdgeColor','none','FontSize',14);

    %plot mean response for each dose with standard error
    errorbar(Doserange,Response,Error,'o','MarkerFaceColor',[0 0 0],'LineWidth',1.5,'Color',[0 0 0],'MarkerSize',9,'HandleVisibility','off');
    hold off;

    %adjust figure
    titletext = append(cellline{c},' ',drug{d});
    title(titletext,'FontSize',22,'FontWeight','normal');
    xlabel('Dose (µM)','FontSize',22);
    ylabel(ylabeltext,'FontSize',22);
    ax = gca;
    %     yticks(0:0.25:1.1); ylim([0 1.2]);
    set(ax,'XLimitMethod','padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
        'XGrid','off','Box','on','Color','none','FontSize',22);

    xticks(xaxislabels);
    xticklabels(string(xaxislabels));

    %save figure
%     filetext = append('DR_plots/',date,'_DR_',cellline{c},'_',drug{d},'_',channel{a},'_Fig6_EC50_combined.svg');
%     saveas(fig, filetext);

    vars = {'Doserange';'Response';'Error';'excl_Dose';'excl_Resp';'excl_Error';'excluded_indices'};
    clear(vars{:});

end %loop d drug

end %function