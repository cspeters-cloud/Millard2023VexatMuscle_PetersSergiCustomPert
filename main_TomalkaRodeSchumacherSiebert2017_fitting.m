%%
% SPDX-FileCopyrightText: 2025 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
%%
clc;
close all;
clear all;

%%
% Parameters
%%
simConfig.runFitting              = 0;
simConfig.generatePlots           = 1;
simConfig.trials                  = [1,2,3];
simConfig.flag_debugFitting       = 1;


modelConfig.fibrilOption     = 'Fibril'; %Fibril or ''
modelConfig.wlcOption        = ''; %WLC or ''
modelConfig.muscleName       = 'EDL';
modelConfig.experimentName   = 'TRSS2017';

fittingConfig.numberOfBisections = 10;


%%
% Set up the files
%%

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

addpath( genpath(projectFolders.parameters)     );
addpath( genpath(projectFolders.curves)         );
addpath( genpath(projectFolders.experiments)    );
addpath( genpath(projectFolders.simulation)     );
addpath( genpath(projectFolders.models)         );
addpath( genpath(projectFolders.postprocessing) );


%%
% Plot parameters
%%
plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  1,...
                            'numberOfVerticalPlotRows',       4,...
                            'flag_fixedPlotWidth',            1,...
                            'plotWidth',                      6,...
                            'plotHeight',                     6,...
                            'flag_usingOctave',               0);

numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = plotLayoutSettings.plotHeight;
flag_usingOctave              = plotLayoutSettings.flag_usingOctave;

plotHorizMarginCm = 3;
plotVertMarginCm  = 3;

pageWidth   = (plotWidth+plotHorizMarginCm)*numberOfHorizontalPlotColumns...
                +plotHorizMarginCm;

pageHeight  = (plotHeight+plotVertMarginCm)*numberOfVerticalPlotRows...
                +plotVertMarginCm;

plotConfigGeneric;

n=0.75;
lineColorsTRSS2017(1,:) = [0,0,0].*(1-n)            + [1,1,1].*(n);
lineColorsTRSS2017(2,:) = ([44,133,229]./255).*(1-n)+ [1,1,1].*(n);
lineColorsTRSS2017(3,:) = ([63,181,175]./255).*(1-n)+ [1,1,1].*(n);

lineColorsSimTRSS2017(1,:) = [0,0,0];
lineColorsSimTRSS2017(2,:) = [44,133,229]./255;
lineColorsSimTRSS2017(3,:) = [63,181,175]./255;

%%
% Load the models parameters
%%

nTrials = length(simConfig.trials);

ratFibrilModels(nTrials) = struct('musculotendon',[],...
                            'sarcomere',[],...
                            'falData',[],...
                            'fpeData',[],...
                            'curves',[],...
                            'fitting',[]);

for idxTrial = simConfig.trials

    ratFibrilModels(idxTrial) = struct('musculotendon',[],...
                            'sarcomere',[],...
                            'falData',[],...
                            'fpeData',[],...
                            'curves',[],...
                            'fitting',[]);
    
    fileName = [    'rat',modelConfig.experimentName,...
                    modelConfig.muscleName,...
                    modelConfig.fibrilOption,...
                    'ActiveTitin',...
                    modelConfig.wlcOption,...
                    '_',num2str(idxTrial),'.mat'];

    filePathRatMuscle = ...
        fullfile(projectFolders.output_structs_FittedModels,fileName);

    tmpModel   = load(filePathRatMuscle);
    
    modelFields = fields(ratFibrilModels(idxTrial));

    for idxField = 1:1:length(modelFields)
        ratFibrilModels(idxTrial).(modelFields{idxField}) = ...
            tmpModel.ratMuscleModelParameters.(modelFields{idxField});
    end

    %
    % Update the default settings to be consistent with a skinned fiber
    %
    ratFibrilModels(idxTrial).sarcomere.scaleECM            =0;
    ratFibrilModels(idxTrial).sarcomere.scaleTitinProximal  =1;
    ratFibrilModels(idxTrial).sarcomere.scaleTitinDistal    =1;
    
    ratFibrilModels(idxTrial).sarcomere.normLengthTitinActinBondMinimum = 0.;
    normLengthTitinActinBondMinimum = ...
        ratFibrilModels(idxTrial).sarcomere.normLengthTitinActinBondMinimum;
    %fprintf('%1.3e norm-length-titin-actin-bond-minimum\n',...
    %        normLengthTitinActinBondMinimum);
    
    %
    % Use the fast-shortening, slow-lengthening model
    %    
    responseTimeScaling=10;
    
    ratFibrilModels(idxTrial).sarcomere.slidingTimeConstantBlendingParameter = 0.01;
    
    ratFibrilModels(idxTrial).sarcomere.slidingTimeConstantLengthening= ...
        ratFibrilModels(idxTrial).sarcomere.slidingTimeConstant...
        *responseTimeScaling;
    
    ratFibrilModels(idxTrial).sarcomere.slidingTimeConstantShortening= ...
        ratFibrilModels(idxTrial).sarcomere.slidingTimeConstant;
    
    ratFibrilModels(idxTrial).sarcomere.useVariableSlidingTimeConstant = 1;

end

%%
% Reference data
%%
[expData, expIndices] = ...
        loadRatSkeletalMuscleData(projectFolders);

expTRSS2017 = expData(expIndices.index_TRSS2017);
%%
% Fitting the active-force-length relaton
%%

if(simConfig.runFitting==1)
    fidFitting = fopen(fullfile(projectFolders.output_structs_TRSS2017,...
                        'fittingLog.txt'),'w');

    figDebugFitting = figure;

    dfalN = 0.2;

    lceOptMdl   = ratFibrilModels(1).musculotendon.optimalFiberLength;
    vmax        = ratFibrilModels(1).musculotendon.maximumNormalizedFiberVelocity;
    lceOptData  = min(expTRSS2017.activeLengtheningData(3).x);

    expfl.lce = [];
    expfl.fN   = []; 

    mdlfl.lce = [];
    mdlfl.fN   = [];

    assert(ratFibrilModels(1).curves.useCalibratedCurves==1,...
           'Error: the calibrated curves should be used');

    paramsFlN = ratFibrilModels(1).curves.activeForceLengthCurve; 

    for idxTrial = simConfig.trials
        lce  = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
        fN    = expTRSS2017.activeLengtheningData(idxTrial).y(1,1);
        

%         if(abs(lceN-1)<1e-3)
%             fN=1;
%         end
        expfl.lce  = [expfl.lce;lce];
        expfl.fN    = [expfl.fN;fN];   

        fNMdl = calcBezierYFcnXDerivative(lce./lceOptMdl,paramsFlN,0);
        mdlfl.lce  = [mdlfl.lce;lce];
        mdlfl.fN    = [mdlfl.fN;fNMdl];
                
    end

    argBest  = 0.2;
    flag_compensateForCrossbridgeStiffness = 0;
    [errFlN,curveL] = calcErrorTRSS2017ForceLengthRelationAscendingLimb(argBest,...
                        expfl.lce,expfl.fN,...
                        ratFibrilModels(1).sarcomere,...
                        ratFibrilModels(1).musculotendon,...
                        flag_compensateForCrossbridgeStiffness);
    errBest = sqrt(mean(errFlN.^2));

    fprintf('%1.2e\tfitting: fal rmse (start)\n',errBest);
    fprintf('%e\tfal-asc offset (start)\n',argBest);

    fprintf(fidFitting,'%1.2e\tfitting: fal rmse (start)\n',errBest);
    fprintf(fidFitting,'%e\tfal-asc offset (start)\n',argBest);    

    argDelta = argBest/2;

    for i=1:1:fittingConfig.numberOfBisections

        argL = argBest-argDelta;
        [errL,curveL] = calcErrorTRSS2017ForceLengthRelationAscendingLimb(argL,...
                            expfl.lce,expfl.fN,...
                            ratFibrilModels(1).sarcomere,...
                            ratFibrilModels(1).musculotendon,...
                            flag_compensateForCrossbridgeStiffness);
        errLMag = sqrt(mean(errL.^2));

        argR = argBest+argDelta;
        [errR,curveR] = calcErrorTRSS2017ForceLengthRelationAscendingLimb(argR,...
                            expfl.lce,expfl.fN,...
                            ratFibrilModels(1).sarcomere,...
                            ratFibrilModels(1).musculotendon,...
                            flag_compensateForCrossbridgeStiffness);
        errRMag = sqrt(mean(errR.^2));        

        if(errLMag < errBest && errLMag <= errRMag )
            argBest=argL;
            errBest=errLMag;
            curveBest=curveL;
        elseif(errRMag < errBest && errRMag < errLMag)
            argBest=argR;
            errBest=errRMag;
            curveBest=curveR;            
        end

        argDelta=argDelta*0.5;

    end

    fprintf('%1.2e\tfitting: fal rmse (end)\n',errBest);
    fprintf('%e\tfal-asc offset (end)\n\n',argBest);

    fprintf(fidFitting,'%1.2e\tfitting: fal rmse (end)\n\n',errBest);
    fprintf(fidFitting,'%e\tfal-asc offset (end)\n\n',argBest);
       



    %
    % Update all of the models
    %
    for i=1:1:simConfig.trials

        flag_compensateForCrossbridgeStiffness=0;
        [flNError, falCurve] = ...
            calcErrorTRSS2017ForceLengthRelationAscendingLimb(...
                argBest, expfl.lce,expfl.fN,...
                ratFibrilModels(i).sarcomere,...
                ratFibrilModels(i).musculotendon,...
                flag_compensateForCrossbridgeStiffness);

        ratFibrilModels(i).curves.activeForceLengthCurve=falCurve;

        flag_compensateForCrossbridgeStiffness=1;
        [flNErrorCal, falCurveCal] = ...
            calcErrorTRSS2017ForceLengthRelationAscendingLimb(...
                argBest, expfl.lce,expfl.fN,...
                ratFibrilModels(i).sarcomere,...
                ratFibrilModels(i).musculotendon,...
                flag_compensateForCrossbridgeStiffness);

        ratFibrilModels(i).curves.activeForceLengthCalibratedCurve=falCurveCal;
        
    end


    if(simConfig.flag_debugFitting==1)

        subplot('Position',reshape(subPlotPanel(1,1,:),1,4));

        sampleFlN=calcBezierYFcnXCurveSampleVector(curveBest,100,curveBest.xEnd);
        plot(sampleFlN.x.*lceOptMdl,sampleFlN.y,'-',...
             'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M)$$');
        hold on;


        for idxTrial = simConfig.trials
            plot(expfl.lce(idxTrial),...
                 expfl.fN(idxTrial),...
                 'x','Color',lineColorsTRSS2017(idxTrial,:),...
                 'MarkerFaceColor',lineColorsTRSS2017(idxTrial,:),...
                 'DisplayName',...
                 expTRSS2017.activeLengtheningData(idxTrial).seriesName);
            hold on;
        end
        

        xlabel(expTRSS2017.activeLengtheningData(3).xName);
        ylabel(expTRSS2017.activeLengtheningData(3).yName);
        title('Force-Length Relation');  
        box off;
        hold on;

    end    

end


%%
% Fitting the properties of the model for the very first 100 ms of data
%%
if(simConfig.runFitting==1)

    fittingFraction = 1/8;
    npts = round(200*fittingFraction);

    lceOptMdl   = ratFibrilModels(1).musculotendon.optimalFiberLength;
    lceOptData  = min(expTRSS2017.activeLengtheningData(3).x);

    if(simConfig.flag_debugFitting==1)

        for idxTrial=simConfig.trials
            if(simConfig.flag_debugFitting==1)
                figure(figDebugFitting);
                subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
                for idxTrial = simConfig.trials
                    plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                         expTRSS2017.activeLengtheningData(idxTrial).y,...
                         '-','Color',lineColorsTRSS2017(idxTrial,:),...
                         'DisplayName',...
                         expTRSS2017.activeLengtheningData(idxTrial).seriesName);
                    hold on;
                end
                xlabel(expTRSS2017.activeLengtheningData(3).xName);
                ylabel(expTRSS2017.activeLengtheningData(3).yName);
                title(expTRSS2017.activeLengtheningData(3).title);  
                box off;
                hold on;
            end
        end
    end

    %
    % Fit the lengthening time constant
    %
    
    bestValue = 20;
    deltaValue=bestValue*0.5;

    optParams.name = 'responseTimeScaling';
    optParams.value=bestValue;
    simConfig.flag_debugFitting=0;

    [bestError,figDebugFitting,ratFibrilModelsUpd,benchRecord] ...
        = calcErrorTRSS2017RampFraction(optParams,...
                   fittingFraction, npts, ...
                   ratFibrilModels, expTRSS2017,simConfig,...
                   figDebugFitting,subPlotPanel,lineColorsSimTRSS2017);

    fprintf('%1.2e\tfitting: sliding time-constant (start)\n',bestError);
    fprintf('%e\t sliding-time constant scaling (end)\n\n',bestValue);

    fprintf(fidFitting,'%1.2e\tfitting: sliding time-constant (start)\n',bestError);
    fprintf(fidFitting,'%e\t sliding-time constant scaling (start)\n\n',bestValue);



    for i=1:1:fittingConfig.numberOfBisections
        fprintf('%i/%i\n',i,fittingConfig.numberOfBisections);

        optParams.value=bestValue-deltaValue;

        [errorVal,figDebugFitting,ratFibrilModelsUpd,benchRecord] ...
            = calcErrorTRSS2017RampFraction(optParams,...
                   fittingFraction, npts, ...
                   ratFibrilModels, expTRSS2017,simConfig,...
                   figDebugFitting,subPlotPanel,lineColorsSimTRSS2017);

        if(errorVal < bestError)
           bestError=errorVal;  
           bestValue = optParams.value;
           ratFibrilModels=ratFibrilModelsUpd;
        else
            optParams.value=bestValue+deltaValue;
            [errorVal,figDebugFitting,ratFibrilModelsUpd,benchRecord] ...
                = calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModels, expTRSS2017,simConfig,...
                       figDebugFitting,subPlotPanel,lineColorsSimTRSS2017);
            if(errorVal < bestError)
                bestError=errorVal;        
                bestValue = optParams.value;               
                ratFibrilModels=ratFibrilModelsUpd;
            end
        end
        deltaValue=deltaValue*0.5;
    end
    fprintf('%1.2e\tfitting: sliding time-constant (end)\n',bestError);
    fprintf('%e\t sliding-time constant scaling (end)\n\n',bestValue);

    %
    % Scale the stiffness and cross-bridge damping
    %
    bestValue = 1;
    deltaValue=bestValue*0.5;

    optParams.name = 'xeStiffnessDampingScaling';
    optParams.value=bestValue;
    simConfig.flag_debugFitting=0;

    [bestError,figDebugFitting,ratFibrilModelsUpd,benchRecord]...
        = calcErrorTRSS2017RampFraction(optParams,...
                   fittingFraction, npts, ...
                   ratFibrilModels, expTRSS2017,simConfig,...
                   figDebugFitting,subPlotPanel,lineColorsSimTRSS2017);

    fprintf('%1.2e\tfitting: impedance scaling (start)\n',bestError);
    fprintf('%e\t impedance scaling (start)\n\n',bestValue);
    fprintf(fidFitting,'%1.2e\tfitting: impedance scaling (start)\n',bestError);
    fprintf(fidFitting,'%e\t impedance scaling (start)\n\n',bestValue);
    
    for i=1:1:fittingConfig.numberOfBisections
        fprintf('%i/%i\n',i,fittingConfig.numberOfBisections);

        optParams.value=bestValue-deltaValue;

        [errorVal,figDebugFitting,ratFibrilModelsUpd,benchRecord] ...
            = calcErrorTRSS2017RampFraction(optParams,...
                   fittingFraction, npts, ...
                   ratFibrilModels, expTRSS2017,simConfig,...
                   figDebugFitting,subPlotPanel,lineColorsSimTRSS2017);

        if(errorVal < bestError)
           bestError=errorVal;           
           bestValue = optParams.value;           
           ratFibrilModels=ratFibrilModelsUpd;
        else
            optParams.value=bestValue+deltaValue;
            [errorVal,figDebugFitting,ratFibrilModelsUpd,benchRecord] ...
                = calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModels, expTRSS2017,simConfig,...
                       figDebugFitting,subPlotPanel,lineColorsSimTRSS2017);
            if(errorVal < bestError)
                bestError=errorVal;      
                bestValue = optParams.value;               
                ratFibrilModels=ratFibrilModelsUpd;
            end
        end
        deltaValue=deltaValue*0.5;
    end
    fprintf('%1.2e\tfitting: impedance scaling (end)\n',bestError);
    fprintf('%e\t impedance scaling (end)\n\n',bestValue);
    fprintf(fidFitting,'%1.2e\tfitting: impedance scaling (end)\n',bestError);    
    fprintf(fidFitting,'%e\t impedance scaling (end)\n\n',bestValue);
    

    save(fullfile(projectFolders.output_structs_TRSS2017,...
            'ratFibrilModelsFitted.mat'),...
            'ratFibrilModels','-mat');

    %
    % Plot the fitted results
    %    
    fullRampFraction=1;
    npts = round(200*fullRampFraction);
    simConfig.flag_debugFitting=1;

    [errorVal,figDebugFitting,ratFibrilModelsUpd,benchRecordFitted] ...
        = calcErrorTRSS2017RampFraction(optParams,...
               fullRampFraction, npts, ...
               ratFibrilModels, expTRSS2017,simConfig,...
               figDebugFitting,subPlotPanel,lineColorsSimTRSS2017);

    save(fullfile(projectFolders.output_structs_TRSS2017,...
            'benchRecordVexat_TRSS2017_fitted.mat'),...
            'benchRecordFitted','-mat');

    fclose(fidFitting);    

end

if(simConfig.generatePlots==1)
    load(fullfile(projectFolders.output_structs_TRSS2017,...
                    'benchRecordVexat_TRSS2017_fitted.mat'));

    figPub=figure;

    nTrials = length(simConfig.trials);
    titinAnalysis(nTrials)=struct('lce',[],'fceN',[],'fxeN',[],'f2N',[],'k2N',[]);
    for idx=simConfig.trials
        titinAnalysis(idx)= struct('lce',[],'fceN',[],'fxeN',[],'f2N',[],'k2N',[]);
    end


    expfl.lce = [];
    expfl.fN   = []; 

    mdlfl.lce = [];
    mdlfl.fN   = [];

    assert(ratFibrilModels(1).curves.useCalibratedCurves==1,...
           'Error: the calibrated curves should be used');

    paramsFlN = ratFibrilModels(1).curves.activeForceLengthCurve; 
    lceOptMdl = ratFibrilModels(1).musculotendon.optimalFiberLength;    

    for idxTrial = simConfig.trials
        lce  = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
        fN   = expTRSS2017.activeLengtheningData(idxTrial).y(1,1);

        expfl.lce  = [expfl.lce;lce];
        expfl.fN    = [expfl.fN;fN];   

        fNMdl = calcBezierYFcnXDerivative(lce./lceOptMdl,paramsFlN,0);
        mdlfl.lce  = [mdlfl.lce;lce];
        mdlfl.fN    = [mdlfl.fN;fNMdl];
                
    end    

    fceNMax   = 2.6;
    kceNMax   = 6;
    idxKceMin = 12;
    lceNMin   =0.6;
    lceNMax = 1.45;

    subplot('Position',reshape(subPlotPanel(1,1,:),1,4));

        curveBest = ratFibrilModels(idx).curves.activeForceLengthCurve;
        
        sampleFlN=calcBezierYFcnXCurveSampleVector(curveBest,100,curveBest.xEnd);
        plot(sampleFlN.x,sampleFlN.y,'-',...
             'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M)$$');
        hold on;
    
        for idx=simConfig.trials        
            plot(expfl.lce(idx)./lceOptMdl,...
                 expfl.fN(idx),...
                 'o','Color',lineColorsTRSS2017(idx,:),...
                 'MarkerFaceColor',lineColorsTRSS2017(idx,:),...
                 'DisplayName',...
                 expTRSS2017.activeLengtheningData(idx).seriesName);
            hold on;
        end
        xlabel('Norm. Length ($$\ell/\ell_o$$)');
        ylabel(expTRSS2017.activeLengtheningData(3).yName);
        title('Force-Length Relation');  
        xlim([lceNMin,lceNMax]);
        box off;
        hold on;   

    subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
        fill([min(sampleFlN.x);...
              max(sampleFlN.x);...
              fliplr(sampleFlN.x')'],...
              [0;...
               0;...
              fliplr(sampleFlN.y')'],...
              [1,1,1].*0.85,'EdgeAlpha',0);
        hold on;  

    for idx=simConfig.trials
        plot(expTRSS2017.activeLengtheningData(idx).x./lceOptMdl,...
             expTRSS2017.activeLengtheningData(idx).y,...
             '-','Color',lineColorsTRSS2017(idx,:),...
             'DisplayName',...
             expTRSS2017.activeLengtheningData(idx).seriesName,...
             'LineWidth',2);
        hold on;  
    end

    for idx=simConfig.trials     
        subplot('Position',reshape(subPlotPanel(2,1,:),1,4)); 



            txtName = expTRSS2017.activeLengtheningData(idx).seriesName;
            i0=strfind(txtName,'Exp.');
            txtName(1,i0:4)='Sim.';

            plot(benchRecordFitted.normFiberLength(:,idx),...
                 benchRecordFitted.normCrossBridgeForceAlongTendon(:,idx),...
                 '-','Color',lineColorsSimTRSS2017(idx,:),...
                 'DisplayName',...
                 txtName);
            hold on;    
            xlim([lceNMin,lceNMax]);

            box off;
            %xlabel(expTRSS2017.activeLengtheningData(3).xName);
            %ylabel(expTRSS2017.activeLengtheningData(3).yName);
            %ylim([0,fceNMax]);
            %title('Estimated total cross-bridge forces');

        subplot('Position',reshape(subPlotPanel(2,1,:),1,4));                


            titinAnalysis(idx).lce = ...
                expTRSS2017.activeLengtheningData(idx).x;
            titinAnalysis(idx).fceN = ...
                expTRSS2017.activeLengtheningData(idx).y;

            [lceU,iU]=unique(benchRecordFitted.normFiberLength(:,idx).*lceOptMdl);
            fxeNU=benchRecordFitted.normCrossBridgeForceAlongTendon(iU,idx);

            titinAnalysis(idx).fxeN = interp1(lceU,fxeNU,titinAnalysis(idx).lce);
            titinAnalysis(idx).f2N  = titinAnalysis(idx).fceN-titinAnalysis(idx).fxeN;
            
            titinAnalysis(idx).k2N = ...
                calcCentralDifferenceDataSeries(...
                    titinAnalysis(idx).lce,...
                    titinAnalysis(idx).f2N);
            titinAnalysis(idx).k2N = titinAnalysis(idx).k2N.*lceOptMdl;

            plot(titinAnalysis(idx).lce./lceOptMdl, ...
                 titinAnalysis(idx).f2N,...
                 '-','Color',lineColorsSimTRSS2017(idx,:));
            hold on;
            box off;
            xlim([lceNMin,lceNMax]);
            ylim([0,fceNMax]);

            xlabel('Norm. Length ($$\ell/\ell_o$$)');
            ylabel(expTRSS2017.activeLengtheningData(3).yName);
            title('Estimated crossbridge and titin forces');  

        subplot('Position',reshape(subPlotPanel(3,1,:),1,4));  
            idxPos = find(titinAnalysis(idx).f2N>0.05 ...
                       & titinAnalysis(idx).k2N>0.05);

            plot(titinAnalysis(idx).lce(idxPos)./lceOptMdl, ...
                 titinAnalysis(idx).k2N(idxPos),...
                 '-','Color',lineColorsSimTRSS2017(idx,:));
            hold on;
            box off;
            xlim([lceNMin,lceNMax]);
            ylim([0,kceNMax]);

            xlabel('Norm. Length ($$\ell/\ell_o$$)');
            ylabel('Norm. Stiffness ($$(f/f_o)/(\ell/\ell_o)$$)');
            title('Estimated titin stiffness');              

    end
    figure(figPub);    
    configPlotExporter;
    filePath = fullfile(projectFolders.output_plots_TRSS2017,...
                        'fig_Sim_TRSS2017_Pub.pdf');
    print('-dpdf', filePath); 
end