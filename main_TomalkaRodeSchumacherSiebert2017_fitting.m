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
simConfig.runFitting              = 1;
simConfig.generatePlots           = 1;
simConfig.trials                  = [1,2,3];
simConfig.defaultTrialId          = 0;
simConfig.useDefaultModel         = 1;
simConfig.flag_debugFitting       = 1;


modelConfig.fibrilOption     = 'Fibril'; %Fibril or ''
modelConfig.wlcOption        = ''; %WLC or ''
modelConfig.muscleName       = 'EDL';
modelConfig.experimentName   = 'TRSS2017';

fittingConfig.fitFl             =1;
fittingConfig.fitFv             =1;
fittingConfig.fitTimeConstant   =1;
fittingConfig.fitKx             =1;
fittingConfig.fitQ              =1;
fittingConfig.fitf1HNPreload    =1;

fittingConfig.numberOfBisections = 10;
fittingConfig.idxFvKey = 3;
fittingConfig.idxFlKey = 1;
fittingConfig.titin.trials = [1,2,3]; 
% [1]    : will fit just trial 1
% [1,2,3]: will fit to trials 1,2,3
fittingConfig.titin.individuallyFit = 0; 
% 0: a set of titin parameters that best fits all trials will be solved
% 1: a set of titin parameters the best fits each individual trial will be solved  
%
fittingConfig.titin.applyToAllTrials = 1;
    %If a fittingConfig.titin.trials is a single trial, then setting
    %this flag to 1 will apply the fitted parameters to all other trials

pubPlotOptions.useSmoothedStiffnessData=1;
pubPlotOptions.plotRawStiffnessData=0;
pubPlotOptions.stiffnessLowerForceBound = 0.05;

pubPlotOptions.fceNMax   = 2.6;
pubPlotOptions.kceNMax   = 6;
pubPlotOptions.lceNMin   =0.6;
pubPlotOptions.lceNMax = 1.45;

%%
% Analysis outputs
%%
nTrials = length(simConfig.trials);
titinAnalysis(nTrials)=struct('lce',[],'fceN',[],'fxeN',[],'f2N',[],'k2N',[]);
for idx=simConfig.trials
    titinAnalysis(idx)= struct('lce',[],'fceN',[],'fxeN',[],'f2N',[],'k2N',[]);
end

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
                            'numberOfVerticalPlotRows',       5,...
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
plotVertMarginCm  = 2.0;

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
lineColors.simXE(1,:) = cs.grey;%cs.red.*(1-n)+[1,1,1].*n;
n=0.25;
lineColors.simXE(2,:) = cs.grey;%cs.red.*(1-n)+[1,1,1].*n; 
n=0.50;
lineColors.simXE(3,:) = cs.grey;%cs.red.*(1-n)+[1,1,1].*n;

n=0;
lineColors.calcTitinF(1,:) = cs.blue.*(1-n)+[1,1,1].*n;
n=0.25;
lineColors.calcTitinF(2,:) = cs.blue.*(1-n)+[1,1,1].*n;
n=0.50;
lineColors.calcTitinF(3,:) = cs.blue.*(1-n)+[1,1,1].*n;

n=0;
lineColors.calcTitinK(1,:) = lineColors.calcTitinF(1,:);
n=0.25;
lineColors.calcTitinK(2,:) = lineColors.calcTitinF(2,:); 
n=0.50;
lineColors.calcTitinK(3,:) = lineColors.calcTitinF(3,:);

n=0;
lineColors.simF(1,:) = cs.red.*(1-n)+[1,1,1].*n;
n=0.25;
lineColors.simF(2,:) = cs.red.*(1-n)+[1,1,1].*n; 
n=0.50;
lineColors.simF(3,:) = cs.red.*(1-n)+[1,1,1].*n;

n=0;
lineColors.simTitinF(1,:) = cs.magenta.*(1-n)+[1,1,1].*n;
n=0.25;
lineColors.simTitinF(2,:) = cs.magenta.*(1-n)+[1,1,1].*n; 
n=0.50;
lineColors.simTitinF(3,:) = cs.magenta.*(1-n)+[1,1,1].*n;

n=0;
lineColors.simTitinK(1,:) = lineColors.simTitinF(1,:);
n=0.25;
lineColors.simTitinK(2,:) = lineColors.simTitinF(2,:); 
n=0.50;
lineColors.simTitinK(3,:) = lineColors.simTitinF(3,:);

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
    
    if(simConfig.useDefaultModel==1)
        fileName = [    'rat',modelConfig.experimentName,...
                        modelConfig.muscleName,...
                        modelConfig.fibrilOption,...
                        'ActiveTitin',...
                        modelConfig.wlcOption,...
                        '_',num2str(simConfig.defaultTrialId),'.mat'];        
    else
        fileName = [    'rat',modelConfig.experimentName,...
                        modelConfig.muscleName,...
                        modelConfig.fibrilOption,...
                        'ActiveTitin',...
                        modelConfig.wlcOption,...
                        '_',num2str(idxTrial),'.mat'];
    end
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


ratFibrilModelsFitted = ratFibrilModels;
fidFitting = [];

fittingTrialsStr = '';
for i=1:1:length(fittingConfig.titin.trials)
    fittingTrialsStr = [fittingTrialsStr,num2str(fittingConfig.titin.trials(1,i))];
end
if(fittingConfig.titin.individuallyFit==1)
    fittingTrialsStr = [fittingTrialsStr,'i'];
end

if(simConfig.runFitting==1)
    fidFitting = fopen(fullfile(projectFolders.output_structs_TRSS2017,...
                    ['fittingLog_',fittingTrialsStr,'.txt']),'w');
end

%%
% Fitting the active-force-length relaton
%%

if(simConfig.runFitting==1 && fittingConfig.fitFl==1)

    %%
    % fal fitting
    %%

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
        [optError, falCurve] = ...
            calcErrorTRSS2017ForceLengthRelationAscendingLimb(...
                argBest, expfl.lce,expfl.fN,...
                ratFibrilModelsFitted(1).sarcomere,...
                ratFibrilModelsFitted(1).musculotendon,...
                flag_compensateForCrossbridgeStiffness);

        ratFibrilModelsFitted(i).curves.activeForceLengthCurve=falCurve;

        flag_compensateForCrossbridgeStiffness=1;
        [optErrorCal, falCurveCal] = ...
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

if(simConfig.runFitting==1 && fittingConfig.fitFv==1)
    %%
    % fv fitting
    %%

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
        idxKey = expTRSS2017.activeLengtheningData(idxTrial).keyIndices(...
                    fittingConfig.idxFvKey); 
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

if(simConfig.runFitting==1 && fittingConfig.fitTimeConstant==1)
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

    
end


if(simConfig.runFitting==1 && fittingConfig.fitKx==1)
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
    

end

%%
% Fit the parameters of titin to get the closest match possible during
% the trial. Here we will fit two parameters:
%
% Q : the point within the PEVK segment that attaches to actin. A value of
%     0 corresponds to the N2A epitope (the most proximal end) while a 
%     value 1 corresponds to the most distal end of the PEVK segment
%
% f1HNPreload: the preload force that the proximal segment develops
%
% Unlike the previous fitting routines, here we need to:
% 1. Fit to one of the trials, and then apply the parameters to all
% 2. Individually fit each trial
%%

if(simConfig.runFitting==1)
    optParams.exp(3) = struct('name','','value',0,'x',[],'y',[],'dydx',[],...
                              'xLine',[],'yLine',[]);
    
    for idxTrial=simConfig.trials
        idx1 = length(expTRSS2017.activeLengtheningData(idxTrial).x);

        x0 = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
        x1 = expTRSS2017.activeLengtheningData(idxTrial).x(end);
        xStart= x0 + 0.2*(x1-x0);
        idx0 = find(expTRSS2017.activeLengtheningData(idxTrial).x > xStart, 1 );

        expXTmp = expTRSS2017.activeLengtheningData(idxTrial).x(idx0:idx1);
        expYTmp = expTRSS2017.activeLengtheningData(idxTrial).y(idx0:idx1);

        A = [expXTmp, ones(size(expXTmp))];
        b = expYTmp;
        p = (A'*A)\(A'*b);
        expSlope = p(1,1);

        xLine = expXTmp;
        yLine = [xLine, ones(size(xLine))]*p;

        lopt=ratFibrilModelsFitted(idxTrial).musculotendon.optimalFiberLength;
        fiso=ratFibrilModelsFitted(idxTrial).musculotendon.fiso;

        optParams.exp(idxTrial).x = expXTmp;
        optParams.exp(idxTrial).y = expYTmp;
        optParams.exp(idxTrial).dydx = expSlope;
        optParams.exp(idxTrial).xLine=xLine;
        optParams.exp(idxTrial).yLine=yLine;

    end

end

if(simConfig.runFitting==1 && fittingConfig.fitQ==1)
    %
    % Solve for Q to best fit the average slope of the force development
    %
    flagDebug=1;
    simConfigTmp = simConfig;
    figDebugFittingQ=figure;
    
    loops = length(fittingConfig.titin.trials);
    if(fittingConfig.titin.individuallyFit==0)
        loops=1;
    end
    
    Q = 0.5;
    QDelta = 0.5*Q;
    
    optParams.name = 'Q';
    optParams.value=  Q;

    for idxLoop = 1:1:loops
        idxTrial = nan;
        if(fittingConfig.titin.individuallyFit==1)
            idxTrial = simConfig.trials(1,idxLoop);
            simConfigTmp.trials = idxTrial;
        end
        
        simConfigTmp.flag_debugFitting=0;
        fittingFraction=1;
        npts=100;
        [optError,figDebugFittingQ,ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                       figDebugFittingQ,subPlotPanel,lineColors.simTitinK);

        optErrorBest=optError;
        QBest=Q;
        dirMap = [1,-1];

        if(fittingConfig.titin.individuallyFit==1)
            fprintf('%1.2e\tfitting trial %i: Q rmse terminal slope (start)\n',...
                    optErrorBest,idxTrial);
            fprintf('%e\tQ (start)\n',QBest);    
    
            fprintf(fidFitting,...
                '%1.2e\tfitting trial %i: Q rmse terminal slope (start)\n',...
                optErrorBest,idxTrial);
            fprintf(fidFitting,...
                '%e\tQ (start)\n',QBest);    
        else
            fprintf('%1.2e\tfitting all: Q rmse terminal slope (start)\n',...
                    optErrorBest);
            fprintf('%e\tQ (start)\n',QBest);    
    
            fprintf(fidFitting,...
                '%1.2e\tfitting all: Q rmse terminal slope (start)\n',...
                optErrorBest);
            fprintf(fidFitting,...
                '%e\tQ (start)\n',QBest);    
        end

        for i=1:1:fittingConfig.numberOfBisections
            fprintf('%i/%i\n',i,fittingConfig.numberOfBisections);
        
            for j=1:1:length(dirMap)
                Q = QBest + QDelta*dirMap(1,j);
                optParams.value=Q;
                [optError,figDebugFittingQ,ratFibrilModelsFittedUpd,benchRecord] =...
                    calcErrorTRSS2017RampFraction(optParams,...
                               fittingFraction, npts, ...
                               ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                               figDebugFittingQ,subPlotPanel,lineColors.simTitinK);
                if(optError<optErrorBest)
                    optErrorBest=optError;
                    QBest=Q;
                    break;
                end
            end
            
            QDelta=QDelta*0.5;

        end

        fprintf('%1.2e\tfitting: Q rmse terminal slope (end)\n',optErrorBest);
        fprintf('%e\tQ (end)\n',QBest);    
        fprintf(fidFitting,'%1.2e\tfitting: Q rmse terminal slope (end)\n',optErrorBest);
        fprintf(fidFitting,'%e\tQ (end)\n',QBest);    

        %
        % Update the parameter struct
        % 
        optParams.value=QBest;
        [optError,figDebugFittingQ,ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                       figDebugFittingQ,subPlotPanel,lineColors.simTitinK);

        if(fittingConfig.titin.individuallyFit==1)
            ratFibrilModelsFitted(idxTrial)=ratFibrilModelsFittedUpd(idxTrial);
        else
            ratFibrilModelsFitted=ratFibrilModelsFittedUpd;
        end
        
    end

    if(flagDebug==1)
        if(exist('figExpSlope')==1)
            clf(figExpSlope);
        else
            figExpSlope = figure;
        end
        for idxTrial=simConfig.trials
            plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                 expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
            hold on;
            plot(optParams.exp(idxTrial).x,optParams.exp(idxTrial).y,'xk');
            hold on;
            plot(optParams.exp(idxTrial).xLine,optParams.exp(idxTrial).yLine,'-b');
            hold on;
            plot(benchRecord.normFiberLength(:,idxTrial).*lopt,...
                 benchRecord.normFiberForce(:,idxTrial),'-r' );
        end
        xlabel('X');
        ylabel('Y');
        title(sprintf('Trial %i',idxTrial));
    end
    
    %
    % If we are fitting just one of the trials, then update the others
    % to have the same fitted Q value
    %
    if(length(fittingConfig.titin.trials)==1 ...
            && fittingConfig.titin.applyToAllTrials==1)
        for i=1:1:length(ratFibrilModelsFitted)
            if(i ~= fittingConfig.titin.trials(1,1))
                ratFibrilModelsFitted(i)=ratFibrilModelsFitted(fittingConfig.titin.trials);
            end
        end
    end
   
end
%
% Solve for f1HNPreload: the preload that the proximal segment of titin
%                        develops to best match the observed
%                        active-lengthening force profile.
%
if(simConfig.runFitting==1 && fittingConfig.fitf1HNPreload == 1)
    
    flagDebug=1;
    simConfigTmp=simConfig;
    figDebugFittingF1HNPreload=figure;

    loops = length(fittingConfig.titin.trials);
    if(fittingConfig.titin.individuallyFit==0)
        loops=1;
    end


    for idxLoop=1:1:loops

        idxTrial=nan;
        if(fittingConfig.titin.individuallyFit==1)
            idxTrial = simConfig.trials(1,idxLoop);
            simConfigTmp.trials = idxTrial;
        end


        f1HNPreload=0.1;
        f1HNPreloadDelta=f1HNPreload*0.5;

        optParams.name  = 'f1HNPreload';
        optParams.value = f1HNPreload;

        simConfigTmp.flag_debugFitting=0;
        fittingFraction=1;
        npts=100;
        [optError,figDebugFittingF1HNPreload,ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,fittingFraction,npts,...
                       ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                       figDebugFittingF1HNPreload,subPlotPanel,...
                       lineColors.simTitinK);

        optErrorBest=optError;
        f1HNPreloadBest=f1HNPreload;
        dirMap = [1,-1];

        if(fittingConfig.titin.individuallyFit==1)
            fprintf('%1.2e\tfitting trial %i: f1HNPreload rmse (start)\n',...
                    optErrorBest,idxTrial);
            fprintf('%e\tf1HNPreload (start)\n',f1HNPreloadBest);    
            fprintf(fidFitting,...
                    '%1.2e\tfitting trial %i: f1HNPreload rmse (start)\n',...
                    optErrorBest,idxTrial);
            fprintf(fidFitting,...
                    '%e\tf1HNPreload (start)\n',f1HNPreloadBest);    
        else
            fprintf('%1.2e\tfitting all: f1HNPreload rmse (start)\n',...
                    optErrorBest);
            fprintf('%e\tf1HNPreload (start)\n',f1HNPreloadBest);    
            fprintf(fidFitting,...
                    '%1.2e\tfitting all: f1HNPreload rmse (start)\n',...
                    optErrorBest);
            fprintf(fidFitting,...
                    '%e\tf1HNPreload (start)\n',f1HNPreloadBest);    

        end
        for i=1:1:fittingConfig.numberOfBisections
            fprintf('%i/%i\n',i,fittingConfig.numberOfBisections);
        
            for j=1:1:length(dirMap)
                f1HNPreload = f1HNPreloadBest + f1HNPreloadDelta*dirMap(1,j);
                optParams.value=f1HNPreload;
                [optError,figDebugFittingQ,ratFibrilModelsFittedUpd,benchRecord] =...
                    calcErrorTRSS2017RampFraction(optParams,...
                               fittingFraction, npts, ...
                               ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                               figDebugFittingQ,subPlotPanel,lineColors.simTitinK);
                if(optError<optErrorBest)
                    optErrorBest=optError;
                    f1HNPreloadBest=f1HNPreload;
                    break;
                end
            end
            
            f1HNPreloadDelta=f1HNPreloadDelta*0.5;

        end
        

        fprintf(fidFitting,...
            '%1.2e\tfitting trial %i: f1HNPreload rmse (end)\n',...
            optErrorBest,idxTrial);
        fprintf(fidFitting,...
            '%e\tf1HNPreload (end)\n',f1HNPreloadBest);    

        %
        % Update the parameter struct
        % 
        optParams.value=f1HNPreloadBest;
        [optError,figDebugFittingF1HNPreload,ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                       figDebugFittingF1HNPreload,subPlotPanel,lineColors.simTitinK);


        if(fittingConfig.titin.individuallyFit==1)
            ratFibrilModelsFitted(idxTrial)=ratFibrilModelsFittedUpd(idxTrial);
        else
            ratFibrilModelsFitted=ratFibrilModelsFittedUpd;
        end


    end

    if(flagDebug==1)
        if(exist('figExpSlope')==1)
            clf(figExpSlope);
        else
            figExpSlope = figure;
        end
        for idxTrial=simConfig.trials
            plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                 expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
            hold on;
            plot(optParams.exp(idxTrial).x,optParams.exp(idxTrial).y,'xk');
            hold on;
            plot(benchRecord.normFiberLength(:,idxTrial).*lopt,...
                 benchRecord.normFiberForce(:,idxTrial),'-r' );
            plot(optParams.exp(idxTrial).xLine,optParams.exp(idxTrial).yLine,'-b');
            hold on;
        end
        xlabel('X');
        ylabel('Y');
        title(sprintf('Trial %i',idxTrial));
    end

    %
    % If we are fitting just one of the trials, then update the others
    % to have the same fitted Q value
    %
    if(length(fittingConfig.titin.trials)==1 ...
            && fittingConfig.titin.applyToAllTrials==1)
        for i=1:1:length(ratFibrilModelsFitted)
            if(i ~= fittingConfig.titin.trials(1,1))
                ratFibrilModelsFitted(i)=ratFibrilModelsFitted(fittingConfig.titin.trials);
            end
        end
    end
    
end

if(simConfig.runFitting==1)
    fclose(fidFitting);
    
    %
    % Simulate all of the trials with the updated parameters and save the
    % results to file
    %
    optParams.name  = 'simulate';
    optParams.value = nan;
    optParams.expX  = [];
    optParams.expY  = [];
    optParams.expDYDX = nan;

    simConfigTmp=simConfig;
    simConfigTmp.trials =[1,2,3];
    simConfigTmp.flag_debugFitting=0;
    fittingFraction=1;
    npts=100;
    [optError,figDebugFittingF1HNPreload,ratFibrilModelsFittedUpd,benchRecordFitted] =...
        calcErrorTRSS2017RampFraction(optParams,fittingFraction,npts,...
                   ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                   figDebugFittingF1HNPreload,subPlotPanel,...
                   lineColors.simTitinK);

    save(fullfile(projectFolders.output_structs_TRSS2017,...
         ['benchRecordVexat_TRSS2017_fitted',fittingTrialsStr,'.mat']),...
         'benchRecordFitted');
    save(fullfile(projectFolders.output_structs_FittedModels,...
        ['ratTRSS2017EDLFibrilActiveTitinFitted',fittingTrialsStr,'.mat']),...
        'ratFibrilModelsFitted');


end


%%
%
%%

if(simConfig.generatePlots==1)
    load(fullfile(projectFolders.output_structs_TRSS2017,...
            ['benchRecordVexat_TRSS2017_fitted',fittingTrialsStr,'.mat']));
    tmp = load(fullfile(projectFolders.output_structs_FittedModels,...
            ['ratTRSS2017EDLFibrilActiveTitinFitted',fittingTrialsStr,'.mat']));
    ratFibrilModelsFitted=tmp.ratFibrilModelsFitted;



    figPub=figure;

    %
    % Plot the force-length relation
    %


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

    %
    % fal-fitting
    %
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
        xlim([pubPlotOptions.lceNMin,pubPlotOptions.lceNMax]);
        box off;
        hold on;   

    %
    % fv-fitting
    %        
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
            idxKey = expTRSS2017.activeLengtheningData(idx).keyIndices(...
                        fittingConfig.idxFvKey);
            plot(expTRSS2017.activeLengtheningData(idx).x(idxKey)./lceOptMdl,...
                 expTRSS2017.activeLengtheningData(idx).y(idxKey),...
                 'x','Color',[0,0,0],...
                 'HandleVisibility','off');
            hold on; 
        end

        xlabel('Norm. Length ($$\ell/\ell_o$$)');
        ylabel(expTRSS2017.activeLengtheningData(3).yName);
        title('Force-Velocity Curve Fitting');  
        xlim([pubPlotOptions.lceNMin,pubPlotOptions.lceNMax]);
        box off;
        hold on;           

    %
    % Experimental analysis plot
    %

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

            idxKey=expTRSS2017.activeLengtheningData(idx).keyIndices(...
                    fittingConfig.idxFvKey);
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
            

            plot(lceNV(idxValid),...
                fxeNV(idxValid),...
                 '--','Color',lineColors.simXE(idx,:),...
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

            idxKeyF0 = expTRSS2017.activeLengtheningData(idx).keyIndices(...
                            fittingConfig.idxFlKey);
            idxKeyF2 = expTRSS2017.activeLengtheningData(idx).keyIndices(...
                             fittingConfig.idxFvKey);

            plot(expTRSS2017.activeLengtheningData(idx).x(idxKeyF0)./lceOptMdl,...
                 expTRSS2017.activeLengtheningData(idx).y(idxKeyF0),...
                 'o','Color',[0,0,0],'HandleVisibility','off');
            hold on;
            plot(expTRSS2017.activeLengtheningData(idx).x(idxKeyF2)./lceOptMdl,...
                 expTRSS2017.activeLengtheningData(idx).y(idxKeyF2),...
                 'x','Color',[0,0,0],'HandleVisibility','off');
            hold on;
            
            xlim([pubPlotOptions.lceNMin,pubPlotOptions.lceNMax]);

            box off;


        subplot('Position',reshape(subPlotPanel(2,1,:),1,4));                

            idxKey=expTRSS2017.activeLengtheningData(idx).keyIndices(...
                    fittingConfig.idxFvKey);
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
            xlim([pubPlotOptions.lceNMin,pubPlotOptions.lceNMax]);
            ylim([0,pubPlotOptions.fceNMax]);

            xlabel('Norm. Length ($$\ell/\ell_o$$)');
            ylabel(expTRSS2017.activeLengtheningData(3).yName);
            title({'Estimated crossbridge and titin forces','during active-lengthening'});

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
            ylim([0,pubPlotOptions.fceNMax]);        
            if(idx==3)
                legend('Location','NorthWest');
            end
        subplot('Position',reshape(subPlotPanel(3,2,:),1,4));
        
            idxKValid = find(titinAnalysis(idx).f2N >...
                        pubPlotOptions.stiffnessLowerForceBound);

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
            ylim([0,pubPlotOptions.kceNMax]);

            xlabel('Norm. Length Change $$(\ell/\ell_o)-\ell_i$$');
            ylabel('Norm. Slope $$(f/f_o)/(\ell/\ell_o)$$');
            title('Estimated titin stiffness-length relation');  
            if(idx==3)
                legend('Location','SouthEast');
            end

    end
end

%
% Generate plots that compare the simulation results to the experiments
%
if(simConfig.generatePlots==1)
    load(fullfile(projectFolders.output_structs_TRSS2017,...
                    'benchRecordVexat_TRSS2017_fitted.mat'));
    tmp = load(fullfile(projectFolders.output_structs_FittedModels,...
                        'ratTRSS2017EDLFibrilActiveTitinFitted.mat'));
    ratFibrilModelsFitted=tmp.ratFibrilModelsFitted;


    curveFitted = ratFibrilModelsFitted(1).curves.activeForceLengthCurve;    
    sampleFlNFitted=calcBezierYFcnXCurveSampleVector(...
                    curveFitted,100,curveFitted.xEnd);

    subplot('Position',reshape(subPlotPanel(4,1,:),1,4));     
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
        
       plot( benchRecordFitted.normFiberLength(:,idx), ...
             benchRecordFitted.normFiberForce(:,idx),...
             '-','Color',lineColors.simF(idx,:),...
             'DisplayName','Sim: $$f^{*}_i$$',...
             'HandleVisibility',hVis);

       hold on;

       plot( benchRecordFitted.normFiberLength(:,idx), ...
             benchRecordFitted.normDistalTitinForce(:,idx),...
             '-','Color',lineColors.simTitinF(idx,:),...
             'DisplayName','Sim: $$f^{T*}_i$$',...
             'HandleVisibility',hVis);

       hold on;

    end


    box off;
    xlim([pubPlotOptions.lceNMin,pubPlotOptions.lceNMax]);
    ylim([0,pubPlotOptions.fceNMax]); 
    xlabel('Norm. Length ($$\ell/\ell_o$$)');
    ylabel(expTRSS2017.activeLengtheningData(3).yName);
    title({'Simulated and estimated crossbridge','and titin forces during active-lengthening'});


    subplot('Position',reshape(subPlotPanel(5,1,:),1,4));     


        for idx=simConfig.trials

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

        end

        for idx=simConfig.trials

           

           plot( benchRecordFitted.normFiberLength(:,idx)...
                  -benchRecordFitted.normFiberLength(1,idx), ...
                 benchRecordFitted.normDistalTitinForce(:,idx),...
                 '-','Color',lineColors.simTitinF(idx,:),...
                 'DisplayName',sprintf('%s%i%s','Sim: $$f^{T*}_',idx,'$$'),...
                 'HandleVisibility','on');

           hold on;            

           lceNSpan = benchRecordFitted.normFiberLength(end,idx)...
                     -benchRecordFitted.normFiberLength(1,idx);
           text(lceNSpan,...
                benchRecordFitted.normDistalTitinForce(end,idx),...
                 sprintf('%1.1f%s',benchRecordFitted.normDistalTitinForce(end,idx),'$$f_o$$'),...
                 'HorizontalAlignment','left',...
                 'VerticalAlignment','top',...
                 'FontSize',8);
           hold on;


        end    

    box off;
    xlim([0,lceNSpan]);
    ylim([0,pubPlotOptions.fceNMax]);

    xlabel('Norm. Length Change $$(\ell/\ell_o)-\ell_i$$');
    ylabel('Norm. Force ($$f/f_o$$)');
    title('Estimated and simulated titin force-length relation');  

    legend('Location','NorthWest');


    subplot('Position',reshape(subPlotPanel(5,2,:),1,4));     


    for idx=simConfig.trials

        idxKValid = find(titinAnalysis(idx).f2N >...
                    pubPlotOptions.stiffnessLowerForceBound);

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
        if(pubPlotOptions.useSmoothedStiffnessData==0)

            seriesColor = lineColors.calcTitinK(idx,:);

            plot(titinAnalysis(idx).lceN(idxKValid)-titinAnalysis(idx).lceN(1,1), ...
                 titinAnalysis(idx).k2N(idxKValid),...
                 '-','Color',seriesColor,...
                 'DisplayName',...
                 sprintf('%s: %1.1f %s','$$\ell_i$$',...
                        titinAnalysis(idx).lceN(1,1),'$$\ell_o$$'));
            hold on;
            text(lceNSpan,titinAnalysis(idx).k2N(end),...
                 sprintf('%1.1f%s',titinAnalysis(idx).k2N(end),'$$f_o$$'),...
                 'HorizontalAlignment','right',...
                 'VerticalAlignment','bottom',...
                 'FontSize',8);
            hold on;  
         
        end
    end

    for idx=simConfig.trials

       k2NSim = calcCentralDifferenceDataSeries(...
                    benchRecordFitted.normFiberLength(:,idx),...
                    benchRecordFitted.normDistalTitinForce(:,idx));

       k2NSim(isnan(k2NSim)) = 0;
       k2NSim(isinf(k2NSim)) = 0;
       

       plot( benchRecordFitted.normFiberLength(:,idx)...
              -benchRecordFitted.normFiberLength(1,idx), ...
             k2NSim,...
             '-','Color',lineColors.simTitinK(idx,:),...
             'DisplayName',sprintf('%s%i%s','Sim: $$k^{T*}_',idx,'$$'),...
             'HandleVisibility','on');

       hold on;            

       lceNSpan = benchRecordFitted.normFiberLength(end,idx)...
                 -benchRecordFitted.normFiberLength(1,idx);

       text(lceNSpan,...
            k2NSim(end),...
             sprintf('%1.1f%s',k2NSim(end),'$$k_o$$'),...
             'HorizontalAlignment','left',...
             'VerticalAlignment','top',...
             'FontSize',8);
       hold on;

    end    

    box off;
    xlim([0,lceNSpan]);
    ylim([0,pubPlotOptions.kceNMax]);

    xlabel('Norm. Length Change $$(\ell/\ell_o)-\ell_i$$');
    ylabel('Norm. Slope $$(f/f_o)/(\ell/\ell_o)$$');
    title('Estimated titin stiffness-length relation');  
    if(idx==3)
        legend('Location','NorthWest');
    end


end

if(simConfig.generatePlots==1)
    figure(figPub);    
    configPlotExporter;
    filePath = fullfile(projectFolders.output_plots_TRSS2017,...
                        'fig_Sim_TRSS2017_Pub.pdf');
    print('-dpdf', filePath); 
end