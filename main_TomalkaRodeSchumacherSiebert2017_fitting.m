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

pubPlotOptions.useSmoothedStiffnessData=1;
pubPlotOptions.plotRawStiffnessData=0;
pubPlotOptions.stiffnessLowerForceBound = 0.05;
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
plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  2,...
                            'numberOfVerticalPlotRows',       4,...
                            'flag_fixedPlotWidth',            1,...
                            'plotWidth',                      7,...
                            'plotHeight',                     7,...
                            'flag_usingOctave',               0);

numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = plotLayoutSettings.plotHeight;
flag_usingOctave              = plotLayoutSettings.flag_usingOctave;

plotHorizMarginCm = 1.5;
plotVertMarginCm  = 1.5;

pageWidth   = (plotWidth+plotHorizMarginCm)*numberOfHorizontalPlotColumns...
                +plotHorizMarginCm;

pageHeight  = (plotHeight+plotVertMarginCm)*numberOfVerticalPlotRows...
                +plotVertMarginCm;

plotConfigGeneric;

cs = getPaulTolColourSchemes('vibrant');

n=0.75;

lineColors.orig(1,:)=[0,0,0];
lineColors.orig(2,:)=[44,133,229]./255;
lineColors.orig(3,:)=[63,181,175]./255;


n=0.25;
lineColors.exp(1,:) = [0,0,0].*(1-n)+[1,1,1].*n;  
n=0.50;
lineColors.exp(2,:) = [0,0,0].*(1-n)+[1,1,1].*n;           
n=0.75;
lineColors.exp(3,:) = [0,0,0].*(1-n)+[1,1,1].*n; 

n=0;
lineColors.simXE(1,:) = cs.red.*(1-n)+[1,1,1].*n;
n=0.25;
lineColors.simXE(2,:) = cs.red.*(1-n)+[1,1,1].*n; 
n=0.50;
lineColors.simXE(3,:) = cs.red.*(1-n)+[1,1,1].*n;

n=0;
lineColors.calcTitinF(1,:) = cs.blue.*(1-n)+[1,1,1].*n;
n=0.25;
lineColors.calcTitinF(2,:) = cs.blue.*(1-n)+[1,1,1].*n;
n=0.50;
lineColors.calcTitinF(3,:) = cs.blue.*(1-n)+[1,1,1].*n;

n=0;
lineColors.calcTitinK(1,:) = cs.magenta.*(1-n)+[1,1,1].*n;
n=0.33;
lineColors.calcTitinK(2,:) = cs.magenta.*(1-n)+[1,1,1].*n; 
n=0.66;
lineColors.calcTitinK(3,:) = cs.magenta.*(1-n)+[1,1,1].*n;


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


ratFibrilModelsFitted = [];
fidFitting = [];

%%
% Fitting the active-force-length relaton
%%

if(simConfig.runFitting==1)
    ratFibrilModelsFitted=ratFibrilModels;

    fidFitting = fopen(fullfile(projectFolders.output_structs_TRSS2017,...
                        'fittingLog.txt'),'w');

    figDebugFitting = figure;

    dfalN = 0.2;

    lceOptMdl   = ratFibrilModelsFitted(1).musculotendon.optimalFiberLength;
    vmax        = ratFibrilModelsFitted(1).musculotendon.maximumNormalizedFiberVelocity;
    lceOptData  = min(expTRSS2017.activeLengtheningData(3).x);

    expfl.lce = [];
    expfl.fN   = []; 
    
    mdlfl.lce = [];
    mdlfl.fN   = [];

    assert(ratFibrilModelsFitted(1).curves.useCalibratedCurves==1,...
           'Error: the calibrated curves should be used');

    curveFalN   = ratFibrilModelsFitted(1).curves.activeForceLengthCurve; 
    curveFvN    = ratFibrilModelsFitted(1).curves.fiberForceVelocityCurve; 

    for idxTrial = simConfig.trials
        lce  = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
        fN    = expTRSS2017.activeLengtheningData(idxTrial).y(1,1);
        
        expfl.lce  = [expfl.lce;lce];
        expfl.fN    = [expfl.fN;fN];   
        
        fNMdl = calcBezierYFcnXDerivative(lce./lceOptMdl,curveFalN,0);
        mdlfl.lce  = [mdlfl.lce;lce];
        mdlfl.fN    = [mdlfl.fN;fNMdl];                
    end


    argBest  = 0.2;
    flag_compensateForCrossbridgeStiffness = 0;
    [errFlN,curveL] = calcErrorTRSS2017ForceLengthRelationAscendingLimb(argBest,...
                        expfl.lce,expfl.fN,...
                        ratFibrilModelsFitted(1).sarcomere,...
                        ratFibrilModelsFitted(1).musculotendon,...
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
                            ratFibrilModelsFitted(1).sarcomere,...
                            ratFibrilModelsFitted(1).musculotendon,...
                            flag_compensateForCrossbridgeStiffness);
        errLMag = sqrt(mean(errL.^2));

        argR = argBest+argDelta;
        [errR,curveR] = calcErrorTRSS2017ForceLengthRelationAscendingLimb(argR,...
                            expfl.lce,expfl.fN,...
                            ratFibrilModelsFitted(1).sarcomere,...
                            ratFibrilModelsFitted(1).musculotendon,...
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
    for i=simConfig.trials

        flag_compensateForCrossbridgeStiffness=0;
        [flNError, falCurve] = ...
            calcErrorTRSS2017ForceLengthRelationAscendingLimb(...
                argBest, expfl.lce,expfl.fN,...
                ratFibrilModelsFitted(1).sarcomere,...
                ratFibrilModelsFitted(1).musculotendon,...
                flag_compensateForCrossbridgeStiffness);

        ratFibrilModelsFitted(i).curves.activeForceLengthCurve=falCurve;

        flag_compensateForCrossbridgeStiffness=1;
        [flNErrorCal, falCurveCal] = ...
            calcErrorTRSS2017ForceLengthRelationAscendingLimb(...
                argBest, expfl.lce,expfl.fN,...
                ratFibrilModelsFitted(1).sarcomere,...
                ratFibrilModelsFitted(1).musculotendon,...
                flag_compensateForCrossbridgeStiffness);

        ratFibrilModelsFitted(i).curves.activeForceLengthCalibratedCurve=falCurveCal;
        
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
                 'x','Color',lineColors.exp(idxTrial,:),...
                 'MarkerFaceColor',lineColors.exp(idxTrial,:),...
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

if(simConfig.runFitting==1)

    lceOptMdl   = ratFibrilModelsFitted(1).musculotendon.optimalFiberLength;
    vmax        = ratFibrilModelsFitted(1).musculotendon.maximumNormalizedFiberVelocity;
    lceOptData  = min(expTRSS2017.activeLengtheningData(3).x);

    mdlfv.lceN = [];
    mdlfv.vceN = [];  
    mdlfv.fN   = [];
    mdlfv.fvN   = [];

    expfv.lceN = [];
    expfv.vceN = [];  
    expfv.fN   = [];
    expfv.fvN   = [];

    for idxTrial = simConfig.trials
        idxKey = expTRSS2017.activeLengtheningData(idxTrial).keyIndices(3); 
        lce  = expTRSS2017.activeLengtheningData(idxTrial).x(idxKey,1);
        lceN = lce/lceOptMdl;
        fN0  = expTRSS2017.activeLengtheningData(idxTrial).y(1,1);
        fN1  = expTRSS2017.activeLengtheningData(idxTrial).y(idxKey,1);
        fvN = fN1/fN0;
        
        expfv.lceN  = [expfv.lceN;lceN];
        expfv.vceN  = [expfv.vceN;0.11];
        expfv.fN  = [expfv.fN;fN1];
        expfv.fvN  = [expfv.fvN;fvN];   
        
        fvNMdl = calcBezierYFcnXDerivative(0.11,curveFvN,0);
        falNMdl = calcBezierYFcnXDerivative(lce/lceOptMdl,curveFalN,0);

        mdlfv.lceN = [mdlfv.lceN;lceN];
        mdlfv.vceN  = [mdlfv.vceN;0.11];
        mdlfv.fN    = [mdlfv.fN;fvNMdl*falNMdl];
        mdlfv.fvN    = [mdlfv.fvN;fvNMdl];                
    end

    arg = 0;
    argBest=arg;
    delta = 0.1;
    

    [errBest, fvCurveBest] = ...
        calcErrorTRSS2017ForceVelocityRelation(...
            arg, expfv, ...
            ratFibrilModelsFitted(1).curves.activeForceLengthCurve,...
            ratFibrilModelsFitted(1).musculotendon);

    fprintf('%1.2e\tfitting: fv rmse (start)\n',errBest);
    fprintf('%e\tfv ecc offset (start)\n',argBest);

    fprintf(fidFitting,'%1.2e\tfitting: fv rmse (start)\n',errBest);
    fprintf(fidFitting,'%e\tfv ecc offset (start)\n',argBest);


    for i=1:1:fittingConfig.numberOfBisections
        arg = argBest-delta;
        [fvNError, fvCurve] = ...
        calcErrorTRSS2017ForceVelocityRelation(...
            arg, expfv, ...
            ratFibrilModelsFitted(1).curves.activeForceLengthCurve,...
            ratFibrilModelsFitted(1).musculotendon);
        if(fvNError<errBest)
            argBest=arg;
            errBest=fvNError;
            fvCurveBest=fvCurve;
        else
            arg = argBest+delta;
            [fvNError, fvCurve] = ...
            calcErrorTRSS2017ForceVelocityRelation(...
                arg, expfv, ...
                ratFibrilModelsFitted(1).curves.activeForceLengthCurve,...
                ratFibrilModelsFitted(1).musculotendon);
            if(fvNError<errBest)
                argBest=arg;
                errBest=fvNError;
                fvCurveBest=fvCurve;
            end
        end

        delta=delta*0.5;
    end

    fprintf('%1.2e\tfitting: fv rmse (end)\n',errBest);
    fprintf('%e\tfv ecc offset (end)\n',argBest);

    fprintf(fidFitting,'%1.2e\tfitting: fv rmse (end)\n',errBest);
    fprintf(fidFitting,'%e\tfv ecc offset (end)\n',argBest);
    

    %
    % Update all of the models
    %

    for i=simConfig.trials

        ratFibrilModelsFitted(i).curves.fiberForceVelocityCurve=fvCurveBest;

        ratFibrilModelsFitted(i).curves.fiberForceVelocityInverseCurve = ...
            createInverseCurve(fvCurveBest);
        
        fvCalCurve=fvCurveBest;

        fvCalCurve.xpts = fvCalCurve.xpts ...
          .*ratFibrilModelsFitted(i).sarcomere.forceVelocityCalibrationFactor;
        
        fvCalCurve.xEnd = fvCalCurve.xEnd ...
          .*ratFibrilModelsFitted(i).sarcomere.forceVelocityCalibrationFactor;
        
        fvCalCurve.dydxEnd = fvCalCurve.dydxEnd ...
          ./ratFibrilModelsFitted(i).sarcomere.forceVelocityCalibrationFactor;        

        ratFibrilModelsFitted(i).curves.fiberForceVelocityCalibratedCurve=...
            fvCalCurve;

        ratFibrilModelsFitted(i).curves.fiberForceVelocityCalibratedInverseCurve=...
            createInverseCurve(fvCalCurve);

    end   
    here=1;
end



%%
% Fitting the properties of the model for the very first 100 ms of data
%%
if(simConfig.runFitting==1)

    fittingFraction = 1/8;
    npts = round(200*fittingFraction);

    lceOptMdl   = ratFibrilModelsFitted(1).musculotendon.optimalFiberLength;
    lceOptData  = min(expTRSS2017.activeLengtheningData(3).x);

    if(simConfig.flag_debugFitting==1)

        for idxTrial=simConfig.trials
            if(simConfig.flag_debugFitting==1)
                figure(figDebugFitting);
                subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
                for idxTrial = simConfig.trials
                    plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                         expTRSS2017.activeLengtheningData(idxTrial).y,...
                         '-','Color',lineColors.exp(idxTrial,:),...
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
                   ratFibrilModelsFitted, expTRSS2017,simConfig,...
                   figDebugFitting,subPlotPanel,lineColors.simXE);

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
                   ratFibrilModelsFitted, expTRSS2017,simConfig,...
                   figDebugFitting,subPlotPanel,lineColors.simXE);

        if(errorVal < bestError)
           bestError=errorVal;  
           bestValue = optParams.value;
           ratFibrilModelsFitted=ratFibrilModelsUpd;
        else
            optParams.value=bestValue+deltaValue;
            [errorVal,figDebugFitting,ratFibrilModelsUpd,benchRecord] ...
                = calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfig,...
                       figDebugFitting,subPlotPanel,lineColors.simXE);
            if(errorVal < bestError)
                bestError=errorVal;        
                bestValue = optParams.value;               
                ratFibrilModelsFitted=ratFibrilModelsUpd;
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
                   ratFibrilModelsFitted, expTRSS2017,simConfig,...
                   figDebugFitting,subPlotPanel,lineColors.simXE);

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
                   ratFibrilModelsFitted, expTRSS2017,simConfig,...
                   figDebugFitting,subPlotPanel,lineColors.simXE);

        if(errorVal < bestError)
           bestError=errorVal;           
           bestValue = optParams.value;           
           ratFibrilModelsFitted=ratFibrilModelsUpd;
        else
            optParams.value=bestValue+deltaValue;
            [errorVal,figDebugFitting,ratFibrilModelsUpd,benchRecord] ...
                = calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfig,...
                       figDebugFitting,subPlotPanel,lineColors.simXE);
            if(errorVal < bestError)
                bestError=errorVal;      
                bestValue = optParams.value;               
                ratFibrilModelsFitted=ratFibrilModelsUpd;
            end
        end
        deltaValue=deltaValue*0.5;
    end
    fprintf('%1.2e\tfitting: impedance scaling (end)\n',bestError);
    fprintf('%e\t impedance scaling (end)\n\n',bestValue);
    fprintf(fidFitting,'%1.2e\tfitting: impedance scaling (end)\n',bestError);    
    fprintf(fidFitting,'%e\t impedance scaling (end)\n\n',bestValue);
    

    save(fullfile(projectFolders.output_structs_FittedModels,...
            'ratTRSS2017EDLFibrilActiveTitinFitted.mat'),...
            'ratFibrilModelsFitted','-mat');

    %
    % Plot the fitted results
    %    
    fullRampFraction=1;
    npts = round(200*fullRampFraction);
    simConfig.flag_debugFitting=1;

    [errorVal,figDebugFitting,ratFibrilModelsUpd,benchRecordFitted] ...
        = calcErrorTRSS2017RampFraction(optParams,...
               fullRampFraction, npts, ...
               ratFibrilModelsFitted, expTRSS2017,simConfig,...
               figDebugFitting,subPlotPanel,lineColors.simXE);

    save(fullfile(projectFolders.output_structs_TRSS2017,...
            'benchRecordVexat_TRSS2017_fitted.mat'),...
            'benchRecordFitted','-mat');

    fclose(fidFitting);    

end

%%
%
%%

if(simConfig.generatePlots==1)
    load(fullfile(projectFolders.output_structs_TRSS2017,...
                    'benchRecordVexat_TRSS2017_fitted.mat'));
    tmp = load(fullfile(projectFolders.output_structs_FittedModels,...
                        'ratTRSS2017EDLFibrilActiveTitinFitted.mat'));
    ratFibrilModelsFitted=tmp.ratFibrilModelsFitted;



    figPub=figure;

    %
    % Plot the force-length relation
    %
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

    curveFalN = ratFibrilModels(1).curves.activeForceLengthCurve; 
    lceOptMdl = ratFibrilModels(1).musculotendon.optimalFiberLength;    

    for idxTrial = simConfig.trials
        lce  = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
        fN   = expTRSS2017.activeLengtheningData(idxTrial).y(1,1);

        expfl.lce  = [expfl.lce;lce];
        expfl.fN    = [expfl.fN;fN];   

        fNMdl = calcBezierYFcnXDerivative(lce./lceOptMdl,curveFalN,0);
        mdlfl.lce  = [mdlfl.lce;lce];
        mdlfl.fN    = [mdlfl.fN;fNMdl];
                
    end    

    fceNMax   = 2.6;
    kceNMax   = 6;
    idxKceMin = 12;
    lceNMin   =0.6;
    lceNMax = 1.45;

    subplot('Position',reshape(subPlotPanel(1,1,:),1,4));

        curveFitted = ratFibrilModelsFitted(1).curves.activeForceLengthCurve;    
        curveOrig = ratFibrilModels(1).curves.activeForceLengthCurve;

        sampleFlNOrig=calcBezierYFcnXCurveSampleVector(...
                        curveOrig,100,curveOrig.xEnd);
        
        sampleFlNFitted=calcBezierYFcnXCurveSampleVector(...
                        curveFitted,100,curveFitted.xEnd);

        plot(sampleFlNOrig.x,sampleFlNOrig.y,'--',...
             'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M)$$');
        hold on;
        plot(sampleFlNFitted.x,sampleFlNFitted.y,'-',...
             'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M) (fitted)$$');
        hold on;
    
        for idx=simConfig.trials        
            plot(expfl.lce(idx)./lceOptMdl,...
                 expfl.fN(idx),...
                 '+','Color',[0,0,0],...
                 'MarkerFaceColor',lineColors.exp(idx,:),...
                 'DisplayName',...
                 expTRSS2017.activeLengtheningData(idx).seriesName);
            hold on;
        end
        xlabel('Norm. Length ($$\ell/\ell_o$$)');
        ylabel(expTRSS2017.activeLengtheningData(3).yName);
        title('Force-Length Relation Fitting');  
        xlim([lceNMin,lceNMax]);
        box off;
        hold on;   

    subplot('Position',reshape(subPlotPanel(1,2,:),1,4));

        plot(sampleFlNFitted.x,sampleFlNFitted.y,'-',...
             'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M)$$');
        hold on;        

        curveFvN = ratFibrilModels(1).curves.fiberForceVelocityCurve;
        curveFvNFitted = ratFibrilModelsFitted(1).curves.fiberForceVelocityCurve;

        fvN = calcBezierYFcnXDerivative(0.11,curveFvN,0);
        fvNFitted = calcBezierYFcnXDerivative(0.11,curveFvNFitted,0);
        
        plot(sampleFlNFitted.x,sampleFlNFitted.y.*fvN,'--',...
             'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M)f^V(v^M)$$');
        hold on;  
        plot(sampleFlNFitted.x,sampleFlNFitted.y.*fvNFitted,'-',...
             'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M)f^V(v^M)$$ (fitted)');
        hold on;  

        for idx=simConfig.trials
            hdlVis='off';
            if(idx==1)
                hdlVis = 'on';
            end

            plot(expTRSS2017.activeLengtheningData(idx).x./lceOptMdl,...
                 expTRSS2017.activeLengtheningData(idx).y,...
                 '-','Color',lineColors.exp(idx,:),...
                 'DisplayName','$$f^{EXP}$$ TRSS2017',...
                 'HandleVisibility',hdlVis);
            hold on; 
            idxKey = expTRSS2017.activeLengtheningData(idx).keyIndices(1,3);
            plot(expTRSS2017.activeLengtheningData(idx).x(idxKey)./lceOptMdl,...
                 expTRSS2017.activeLengtheningData(idx).y(idxKey),...
                 'x','Color',[0,0,0],...
                 'HandleVisibility','off');
            hold on; 
        end

        xlabel('Norm. Length ($$\ell/\ell_o$$)');
        ylabel(expTRSS2017.activeLengtheningData(3).yName);
        title('Force-Velocity Curve Fitting');  
        xlim([lceNMin,lceNMax]);
        box off;
        hold on;           

    subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
        fill([min(sampleFlNFitted.x);...
              max(sampleFlNFitted.x);...
              fliplr(sampleFlNFitted.x')'],...
              [0;...
               0;...
              fliplr(sampleFlNFitted.y')'],...
              [1,1,1].*0.85,'EdgeAlpha',0,...
              'HandleVisibility','off');
        hold on;  

    for idx=simConfig.trials
        hVis='off';
        if(idx==1)
            hVis='on';
        end
        plot(expTRSS2017.activeLengtheningData(idx).x./lceOptMdl,...
             expTRSS2017.activeLengtheningData(idx).y,...
             '-','Color',lineColors.exp(idx,:),...
             'DisplayName','TRSS2017: $$f^{EXP}_i$$',...
             'HandleVisibility',hVis);
        hold on;  
        text(expTRSS2017.activeLengtheningData(idx).x(end)./lceOptMdl,...
             expTRSS2017.activeLengtheningData(idx).y(end),...
             sprintf('%s%i%s','$$f^{EXP}_{',idx,'}$$'),...
             'HorizontalAlignment','left',...
             'VerticalAlignment','baseline',...
             'FontSize',8);
        hold on;
    end



    for idx=simConfig.trials     
        subplot('Position',reshape(subPlotPanel(2,1,:),1,4)); 

            idxKey=expTRSS2017.activeLengtheningData(idx).keyIndices(3);
            lceNKey = expTRSS2017.activeLengtheningData(idx).x(idxKey)/lceOptMdl;
            fceNKey = expTRSS2017.activeLengtheningData(idx).y(idxKey); 

            fvCurve = ratFibrilModelsFitted(idx).curves.fiberForceVelocityCurve;
            falCurve = ratFibrilModelsFitted(idx).curves.activeForceLengthCurve;
            lceNV = benchRecordFitted.normFiberLength(:,idx);
            vceNV = benchRecordFitted.normFiberVelocity(:,idx);



            fvNV    = zeros(size(lceNV));
            falNV   = zeros(size(lceNV));
            for i=1:1:length(lceNV)
                falNV(i,1)=calcBezierYFcnXDerivative(lceNV(i,1),falCurve,0);
                fvNV(i,1)=calcBezierYFcnXDerivative(vceNV(i,1),fvCurve,0);
            end
            fxeNV = fvNV.*falNV;

            txtName = expTRSS2017.activeLengtheningData(idx).seriesName;
            i0=strfind(txtName,'Exp.');
            txtName(1,i0:4)='Sim.';


            idxValid = find(lceNV >= lceNKey);
    
            hdlVis='off';
            if(idx==1)
                hdlVis='on';                
            end
            
%             fill([lceNV(idxValid(1,1));lceNV(idxValid(end));fliplr(lceNV(idxValid)')'],...
%                  [0;0;fliplr(fxeNV(idxValid)')'],...
%                  lineColors.simXE(idx,:),...
%                  'EdgeColor','none',...
%                  'DisplayName','Crossbridges: $$f^{XE}=f^L(\ell)f^V(v)$$')

            plot(lceNV(idxValid),...
                fxeNV(idxValid),...
                 '-','Color',lineColors.simXE(idx,:),...
                 'DisplayName','Crossbridges: $$f^{XE}_i=f^L(\ell)f^V(v)$$',...
                 'HandleVisibility',hdlVis);
            hold on;   
            text(lceNV(idxValid(1,1)),...
                fxeNV(idxValid(1,1)),...
                sprintf('%s%i%s','$$f^{XE}_{',idx,'}$$'),...
                'FontSize',8,...
                'HorizontalAlignment','left',...
                'VerticalAlignment','top');
            hold on

            idxKeyF0 = expTRSS2017.activeLengtheningData(idx).keyIndices(1);
            idxKeyF2 = expTRSS2017.activeLengtheningData(idx).keyIndices(3);

            plot(expTRSS2017.activeLengtheningData(idx).x(idxKeyF0)./lceOptMdl,...
                 expTRSS2017.activeLengtheningData(idx).y(idxKeyF0),...
                 'o','Color',[0,0,0],'HandleVisibility','off');
            hold on;
            plot(expTRSS2017.activeLengtheningData(idx).x(idxKeyF2)./lceOptMdl,...
                 expTRSS2017.activeLengtheningData(idx).y(idxKeyF2),...
                 'x','Color',[0,0,0],'HandleVisibility','off');
            hold on;
            
            xlim([lceNMin,lceNMax]);

            box off;


        subplot('Position',reshape(subPlotPanel(2,1,:),1,4));                

            idxKey=expTRSS2017.activeLengtheningData(idx).keyIndices(3);
            lceNKey = expTRSS2017.activeLengtheningData(idx).x(idxKey)/lceOptMdl;
            fceNKey = expTRSS2017.activeLengtheningData(idx).y(idxKey);           
                  
            lceNVU= unique(lceNV);

            lceNExp = expTRSS2017.activeLengtheningData(idx).x ./ lceOptMdl;
            [lceNExpU,iu] = unique(lceNExp);
            fceNExp = expTRSS2017.activeLengtheningData(idx).y;
            fceNExpU = fceNExp(iu);

            titinAnalysis(idx).lceN = lceNVU;
            titinAnalysis(idx).fceN = interp1(lceNExpU,fceNExpU,lceNVU,'linear');

            

            [lceNU,iU]=unique(lceNV);
            fxeNU=fxeNV(iU);

            idxZeroTitin = find(titinAnalysis(idx).lceN <= lceNKey);  
            titinAnalysis(idx).fxeN = interp1(lceNU,fxeNU,titinAnalysis(idx).lceN);
            titinAnalysis(idx).f2N  = titinAnalysis(idx).fceN-titinAnalysis(idx).fxeN;
            titinAnalysis(idx).f2N(idxZeroTitin)=0;
            titinAnalysis(idx).f2N(titinAnalysis(idx).f2N<=0)=0;

            hVis='off';
            if(idx==1)
                hVis='on';
            end
            plot(titinAnalysis(idx).lceN, ...
                 titinAnalysis(idx).f2N,...
                 '-','Color',lineColors.calcTitinF(idx,:),...
                 'DisplayName','Titin: $$f^T_i=f^{EXP}_i-f^{XE}_i$$',...
                 'HandleVisibility',hVis);

            text(titinAnalysis(idx).lceN(end),...
                 titinAnalysis(idx).f2N(end),...
                 sprintf('%s%i%s','$$f^{T}_{',idx,'}$$'),...
                 'HorizontalAlignment','left',...
                 'VerticalAlignment','baseline',...
                 'FontSize',8);
                  
            hold on;
            box off;
            xlim([lceNMin,lceNMax]);
            ylim([0,fceNMax]);

            xlabel('Norm. Length ($$\ell/\ell_o$$)');
            ylabel(expTRSS2017.activeLengtheningData(3).yName);
            title('Estimated crossbridge and titin force distribution');  
            if(idx==3)
                legend('Location','NorthWest');
            end

        subplot('Position',reshape(subPlotPanel(2,2,:),1,4));  

            idxKValid = find(titinAnalysis(idx).f2N >...
                        pubPlotOptions.stiffnessLowerForceBound);

            xSpan = max(titinAnalysis(idx).lceN)-min(titinAnalysis(idx).lceN);
            [sp,f2NS,rho] = spaps(...
                titinAnalysis(idx).lceN(idxKValid),...
                titinAnalysis(idx).f2N(idxKValid),1e-6);

            f2NSFit = fnval(sp,titinAnalysis(idx).lceN(idxKValid));
            titinAnalysis(idx).f2NS=titinAnalysis(idx).f2N;
            titinAnalysis(idx).f2NS(idxKValid)=f2NSFit;


            plot(titinAnalysis(idx).lceN,...
                 titinAnalysis(idx).f2N,...
                 '-','Color',lineColors.calcTitinF(idx,:));
            hold on;
            plot(titinAnalysis(idx).lceN,...
                 titinAnalysis(idx).f2NS,...
                 '-','Color',[1,1,1].*0.75);
            hold on;
            box off;
            xlabel('Norm. Length ($$\ell/\ell_o$$)');
            ylabel('Norm. Force ($$f/f_o$$)');
            title('Calc. Titin forces vs. smoothed version');
                

        subplot('Position',reshape(subPlotPanel(3,1,:),1,4));  

            titinAnalysis(idx).k2N = ...
                calcCentralDifferenceDataSeries(...
                    titinAnalysis(idx).lceN,...
                    titinAnalysis(idx).f2N);

            titinAnalysis(idx).k2NS = ...
                calcCentralDifferenceDataSeries(...
                    titinAnalysis(idx).lceN,...
                    titinAnalysis(idx).f2NS);
            lceNSpan = max(titinAnalysis(idx).lceN)...
                      -min(titinAnalysis(idx).lceN);

            hVis='off';
            if(idx==1)
                hVis='on';
            end

            plot(titinAnalysis(idx).lceN-titinAnalysis(idx).lceN(1,1), ...
                 titinAnalysis(idx).f2N,...
                 '-','Color',lineColors.calcTitinF(idx,:),...
                 'DisplayName',...
                 sprintf('%s%i%s','$$f^{T}_{',idx,'}$$'));
            hold on;

            text(lceNSpan,titinAnalysis(idx).f2N(end),...
                 sprintf('%1.1f%s',titinAnalysis(idx).f2N(end),'$$f_o$$'),...
                 'HorizontalAlignment','right',...
                 'VerticalAlignment','bottom',...
                 'FontSize',8);
            hold on;

            box off;

            xlabel('Norm. Length Change $$(\ell/\ell_o)-\ell_i$$');
            ylabel('Norm. Force ($$f/f_o$$)');
            title('Estimated titin force-length relation');  
            
            xlim([0,lceNSpan]);
            ylim([0,fceNMax]);        
            if(idx==3)
                legend('Location','NorthWest');
            end
        subplot('Position',reshape(subPlotPanel(3,2,:),1,4));

            if(pubPlotOptions.useSmoothedStiffnessData==1)
                plot(titinAnalysis(idx).lceN(idxKValid(2:end))-titinAnalysis(idx).lceN(1,1), ...
                     titinAnalysis(idx).k2NS(idxKValid(2:end)),...
                     '-','Color',lineColors.calcTitinK(idx,:),...
                     'DisplayName',...
                     sprintf('%s%i%s','$$k^{T}_{',idx,'}$$'));
                hold on;
                text(lceNSpan,titinAnalysis(idx).k2NS(end),...
                     sprintf('%1.1f',titinAnalysis(idx).k2NS(end)),...
                     'HorizontalAlignment','right',...
                     'VerticalAlignment','bottom',...
                     'FontSize',8);
                hold on;                
            end
            if(pubPlotOptions.plotRawStiffnessData==1 ...
                    || pubPlotOptions.useSmoothedStiffnessData==0)

                seriesColor = lineColors.calcTitinK(idx,:);
                if(pubPlotOptions.plotRawStiffnessData==1)
                    seriesColor = [1,1,1].*0.5;
                end

                plot(titinAnalysis(idx).lceN(idxKValid)-titinAnalysis(idx).lceN(1,1), ...
                     titinAnalysis(idx).k2N(idxKValid),...
                     '-','Color',seriesColor,...
                     'DisplayName',...
                     sprintf('%s: %1.1f %s','$$\ell_i$$',...
                            titinAnalysis(idx).lceN(1,1),'$$\ell_o$$'));
                hold on;
                text(lceNSpan,titinAnalysis(idx).k2N(end),...
                     sprintf('%1.1f',titinAnalysis(idx).k2N(end)),...
                     'HorizontalAlignment','right',...
                     'VerticalAlignment','bottom',...
                     'FontSize',8);
                hold on;  
             

            end
            
            box off;
            xlim([0,lceNSpan]);
            ylim([0,kceNMax]);

            xlabel('Norm. Length Change $$(\ell/\ell_o)-\ell_i$$');
            ylabel('Norm. Slope $$(f/f_o)/(\ell/\ell_o)$$');
            title('Estimated titin stiffness-length relation');  
            if(idx==3)
                legend('Location','SouthEast');
            end

    end
    figure(figPub);    
    configPlotExporter;
    filePath = fullfile(projectFolders.output_plots_TRSS2017,...
                        'fig_Sim_TRSS2017_Pub.pdf');
    print('-dpdf', filePath); 
end