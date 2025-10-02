function [optError,optErrorValues, figDebugFitting,ratFibrilModelsUpd,benchRecord] =...
    calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModels, expTRSS2017,simConfig,...
                       figDebugFitting,subPlotPanel,lineColorsSimTRSS2017)

    ratFibrilModelsUpd=ratFibrilModels;
    errVec = zeros(length(ratFibrilModelsUpd),npts);
    numberOfSimulations = length(simConfig.trials);
    benchRecord             = [];

    optError=0;
    optErrorValues.x = [];
    optErrorValues.y = [];

    for idxTrial = simConfig.trials

        %Update the model parameters
        switch optParams.name

            case 'responseTimeScaling'
                ratFibrilModelsUpd(idxTrial).sarcomere.slidingTimeConstantLengthening= ...
                    ratFibrilModelsUpd(idxTrial).sarcomere.slidingTimeConstant...
                    *optParams.value;     
            case 'xeStiffnessDampingScaling'
                ratFibrilModelsUpd(idxTrial).sarcomere.normCrossBridgeStiffness = ...
                    ratFibrilModelsUpd(idxTrial).sarcomere.normCrossBridgeStiffness ...
                    *optParams.value;
                ratFibrilModelsUpd(idxTrial).sarcomere.normCrossBridgeDamping = ...
                    ratFibrilModelsUpd(idxTrial).sarcomere.normCrossBridgeDamping ...
                    *optParams.value;
            case 'QToF'
                ratFibrilModelsUpd(idxTrial).sarcomere.normPevkToActinAttachmentPoint = ...
                    optParams.value;

              [ratFibrilModelsUpd(idxTrial).curves.forceLengthProximalTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthProximalTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthDistalTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthDistalTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthIgPTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthIgPTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthPevkTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthPevkTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthIgDTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthIgDTitinInverseCurve] ...
                        = createTitinCurves2022( ...
                             ratFibrilModelsUpd(idxTrial).curves.fiberForceLengthCurve,...                                   
                             ratFibrilModelsUpd(idxTrial).curves.forceLengthCurveSettings,...
                             ratFibrilModelsUpd(idxTrial).curves.forceLengthECMHalfCurve,...
                             ratFibrilModelsUpd(idxTrial).sarcomere,...
                             ratFibrilModelsUpd(idxTrial).musculotendon.name,...
                             ratFibrilModelsUpd(idxTrial).curves.useWLCTitinModel,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_createTwoSidedCurves,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_computeCurveIntegrals,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_useElasticIgD,...
                             ratFibrilModelsUpd(idxTrial).sarcomere.titinModelType,...                                   
                             ratFibrilModelsUpd(idxTrial).curves.flag_useOctave);

            case 'QToK'
                ratFibrilModelsUpd(idxTrial).sarcomere.normPevkToActinAttachmentPoint = ...
                    optParams.value;

              [ratFibrilModelsUpd(idxTrial).curves.forceLengthProximalTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthProximalTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthDistalTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthDistalTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthIgPTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthIgPTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthPevkTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthPevkTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthIgDTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthIgDTitinInverseCurve] ...
                        = createTitinCurves2022( ...
                             ratFibrilModelsUpd(idxTrial).curves.fiberForceLengthCurve,...                                   
                             ratFibrilModelsUpd(idxTrial).curves.forceLengthCurveSettings,...
                             ratFibrilModelsUpd(idxTrial).curves.forceLengthECMHalfCurve,...
                             ratFibrilModelsUpd(idxTrial).sarcomere,...
                             ratFibrilModelsUpd(idxTrial).musculotendon.name,...
                             ratFibrilModelsUpd(idxTrial).curves.useWLCTitinModel,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_createTwoSidedCurves,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_computeCurveIntegrals,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_useElasticIgD,...
                             ratFibrilModelsUpd(idxTrial).sarcomere.titinModelType,...                                   
                             ratFibrilModelsUpd(idxTrial).curves.flag_useOctave);
              
            case 'f1HNPreload'
                ratFibrilModelsUpd(idxTrial).sarcomere.f1HNPreload= ...
                    optParams.value;
            case 'simulate'
                %nothing to do here.
            otherwise
                assert(0,'Error: invalid optParams.name');
        end

        lceOptMdl   = ratFibrilModelsUpd(idxTrial).musculotendon.optimalFiberLength;
        vmax        = ratFibrilModelsUpd(idxTrial).musculotendon.maximumNormalizedFiberVelocity;
        lceOptData  = min(expTRSS2017.activeLengtheningData(3).x);

        rampLengthStart  = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
        rampLengthEnd    = expTRSS2017.activeLengtheningData(idxTrial).x(end,1); 
        
        timeStart       = 0;
        timeRampStart   = 0.1;
        rampVelocity    = 0.11*vmax*lceOptData;
        timeRampEnd     = timeRampStart + ...
                          fittingFraction*(...
                            rampLengthEnd-rampLengthStart)/rampVelocity;         
        timeEnd         = timeRampEnd;


        timeSpan = [timeStart,timeEnd];
    
        timeStimulation = timeStart;
        activation      = 1;
    
        excitationFcn = @(argT)calcStepFunction(argT,...
                          timeStimulation-1,...
                          inf,...
                          activation);
        
        pathLengthFcn = @(argT)calcRampStateSharp(argT,...
                                timeRampStart,timeRampEnd,...
                                rampLengthStart,rampVelocity);   

        activationFcn = ...
            @(argU,argA)calcFirstOrderActivationDerivative(argU,argA, ...
                ratFibrilModelsUpd(idxTrial).sarcomere.activationTimeConstant,...
                ratFibrilModelsUpd(idxTrial).sarcomere.deactivationTimeConstant,0);        

        %
        % Bench config
        %
        benchConfig.npts                  = npts;
        benchConfig.relTol                = 1e-6;
        benchConfig.absTol                = 1e-6;
        benchConfig.minActivation         = 0;
        benchConfig.color0                = [0,0,1].*0.5;
        benchConfig.color1                = [0,0,1];
        
        nStates     = 3;
        labelStates = {'$$\dot{\ell}_{a}$$', '$$\ell_{a}$$', '$$\ell_1$$'};       
        
        benchConfig.numberOfMuscleStates  = nStates;
        benchConfig.stateLabels           = labelStates;
        benchConfig.name                  = '';
        benchConfig.initialState          = [];
        benchConfig.initialActivation     = excitationFcn(0);
        benchConfig.pathFcn               = [];
        benchConfig.excitationFcn         = [];
        benchConfig.activationFcn         = activationFcn; 
        benchConfig.tspan                 = timeSpan;  
        
        benchConfig.useFiberDamping  = 1;
        benchConfig.useElasticTendon = 0;
        benchConfig.damping          = 0.1;
        benchConfig.iterMax          = 100;
        benchConfig.tol              = 1e-6;
        
        loopTolerance = min(benchConfig.relTol,benchConfig.absTol)/100;
        
        %
        % Initialization 
        %
        
        modelConfig = struct( ...
          'iterMax'                 , 100             , ...
          'tol'                     , loopTolerance   , ... 
          'tolInit'                 , sqrt(eps)       , ...
          'minActivation'           , 0.0             , ...
          'useElasticTendon'        , 0 , ...
          'initializeState'         , 0                     );  
    
        modelConfig.initializeState =1;
        activationState0    = [0;excitationFcn(0)];
        pathState0          = pathLengthFcn(0);
        muscleState0        = zeros(nStates,1);
        mtInfo = calcMillard2023VexatMuscleInfo(...
                    activationState0,...
                    pathState0,...
                    muscleState0,...
                    ratFibrilModelsUpd(idxTrial).musculotendon,...
                    ratFibrilModelsUpd(idxTrial).sarcomere,...
                    ratFibrilModelsUpd(idxTrial).curves,...
                    modelConfig);

        muscleState0                = mtInfo.state.value;
        modelConfig.initializeState = 0;           
        
        benchConfig.numberOfMuscleStates = length(muscleState0);
        benchConfig.initialState         = muscleState0;
        
        benchConfig.minimumActivation    = 0;
        benchConfig.name                 = 'Vexat';
        benchConfig.eventFcn             = [];            
        
        calcMillard2023VexatMuscleInfoFcn = ...
             @(activationState1,pathState2,mclState3) ...
             calcMillard2023VexatMuscleInfo(    ...
                activationState1,...
                pathState2,...
                mclState3,...
                ratFibrilModelsUpd(idxTrial).musculotendon,...
                ratFibrilModelsUpd(idxTrial).sarcomere,...
                ratFibrilModelsUpd(idxTrial).curves,...
                modelConfig);
        
        %
        % Run the simulation
        %
        idx                     = idxTrial;
        flag_appendEnergetics   = 0;
        flag_useOctave          = 0;
        
        benchConfig.pathFcn               = pathLengthFcn;
        benchConfig.excitationFcn         = excitationFcn;             
        benchRecord = runPrescribedLengthActivationSimulation2025(...
                                   calcMillard2023VexatMuscleInfoFcn,...
                                   [],...
                                   benchConfig,...
                                   benchRecord,...
                                   idx, ...
                                   numberOfSimulations,...
                                   flag_appendEnergetics,...
                                   flag_useOctave);

        %fprintf('%i / %i\n', idx, numberOfSimulations);  

        switch optParams.name
            case 'responseTimeScaling'
                [expLceU,iq] = unique(expTRSS2017.activeLengtheningData(idxTrial).x);
                expfNU = expTRSS2017.activeLengtheningData(idxTrial).y(iq);
        
                errV = zeros(npts,1);
                for k=1:1:npts
                    lceN = benchRecord.normFiberLength(k,idx).*lceOptMdl;
                    fN   = benchRecord.normFiberForce(k,idx);
                    
                    expfN = interp1(expLceU,...
                                    expfNU,...
                                    lceN);
                    errV(k,1) =(expfN-fN);
                end
                if(isempty(optErrorValues))
                    optErrorValues.x = zeros(npts,3);                    
                    optErrorValues.y = zeros(npts,3);
                end
                optErrorValues.x(:,idx) = ...
                    benchRecord.normFiberLength(:,idx).*lceOptMdl;
                optErrorValues.y(:,idx) = errV;
                optError = optError+sqrt(mean(errV.^2));                
            case 'xeStiffnessDampingScaling'
                [expLceU,iq] = unique(expTRSS2017.activeLengtheningData(idxTrial).x);
                expfNU = expTRSS2017.activeLengtheningData(idxTrial).y(iq);
        
                errV = zeros(npts,1);
                for k=1:1:npts
                    lceN = benchRecord.normFiberLength(k,idx).*lceOptMdl;
                    fN   = benchRecord.normFiberForce(k,idx);
                    
                    expfN = interp1(expLceU,...
                                    expfNU,...
                                    lceN);
                    errV(k,1) = expfN-fN;
                end
                if(isempty(optErrorValues))
                    optErrorValues.x = zeros(npts,3);
                    optErrorValues.y = zeros(npts,3);
                end
                optErrorValues.x(:,idx) = ...
                    benchRecord.normFiberLength(:,idx).*lceOptMdl;
                optErrorValues.y(:,idx) = errV;                
                optError = optError+sqrt(mean(errV.^2));
            case 'QToF'
               lopt = ratFibrilModelsUpd(idxTrial).sarcomere.optimalSarcomereLength;

               x0N = optParams.exp(idxTrial).x(1,1)/lopt;
               x1N = optParams.exp(idxTrial).x(end)/lopt;
               idx0 = find(benchRecord.normFiberLength(:,idx)<x0N,1,'last');
               idx1 = size(benchRecord.normFiberLength,1);

               %
               % Sample the simulation record at the same lengths as the
               % experimental data
               %
               mdlX = optParams.exp(idxTrial).x;               
               mdlY = zeros(size(optParams.exp(idxTrial).x));

               for i=1:1:length(mdlX)
                   mdlY(i,1)=interp1(benchRecord.normFiberLength(idx0:idx1,idx).*lopt,...
                                     benchRecord.normFiberForce(idx0:idx1,idx),...
                                     mdlX(i,1));                   
               end

               if(isempty(optErrorValues))
                 optErrorValues.x = zeros(length(mdlX),3);
                 optErrorValues.y = zeros(length(mdlX),3);
               end
               optErrorValues.x(:,idx) = mdlX;
               optErrorValues.y(:,idx) = (mdlY - optParams.exp(idxTrial).y);               
                
               optError = optError+sqrt(mean((mdlY - optParams.exp(idxTrial).y).^2));                       
            case 'QToK'
               lopt = ratFibrilModelsUpd(idxTrial).sarcomere.optimalSarcomereLength;

               x0N = optParams.exp(idxTrial).x(1,1)/lopt;
               x1N = optParams.exp(idxTrial).x(end)/lopt;
               idx0 = find(benchRecord.normFiberLength(:,idx)<x0N,1,'last');
               idx1 = size(benchRecord.normFiberLength,1);

               %
               % Sample the simulation record at the same lengths as the
               % experimental data
               %
               mdlX = optParams.exp(idxTrial).x;               
               mdlY = zeros(size(optParams.exp(idxTrial).x));

                
               %
               % Fit a line to this data
               %
               A = [mdlX,ones(size(mdlX))];
               b= mdlY;
               p = (A'*A)\(A'*b);
               mdlSlope = p(1,1);

               %
               % Calculate the error as the squared difference between the 
               % two slopes
               %              
               if(isempty(optErrorValues))
                 optErrorValues.x = zeros(1,3);
                 optErrorValues.y = zeros(1,3);
               end
               optErrorValues.x(:,idx) = [min(mdlX);max(mdlX)];                              
               optErrorValues.y(:,idx) = [1;1].*(mdlSlope-optParams.exp(idxTrial).dydx);               

               optError = optError+(mdlSlope-optParams.exp(idxTrial).dydx).^2;
               
            case 'f1HNPreload'
               lopt = ratFibrilModelsUpd(idxTrial).sarcomere.optimalSarcomereLength;

               x0N = optParams.exp(idxTrial).x(1,1)/lopt;
               x1N = optParams.exp(idxTrial).x(end,1)/lopt;
               idx0 = find(benchRecord.normFiberLength(:,idx)<x0N,1,'last');
               idx1 = size(benchRecord.normFiberLength,1);

               %
               % Sample the simulation record at the same lengths as the
               % experimental data
               %
               mdlX = optParams.exp(idxTrial).x;               
               mdlY = zeros(size(optParams.exp(idxTrial).x));

               for i=1:1:length(mdlX)
                   mdlY(i,1)=interp1(benchRecord.normFiberLength(idx0:idx1,idx).*lopt,...
                                     benchRecord.normFiberForce(idx0:idx1,idx),...
                                     mdlX(i,1));                   
               end

               if(isempty(optErrorValues))
                 optErrorValues.x = zeros(length(mdlX),3);
                 optErrorValues.y = zeros(length(mdlX),3);
               end
               optErrorValues.x(:,idx) = mdlX;
               optErrorValues.y(:,idx) = mdlY - optParams.exp(idxTrial).y;               

               optError = optError+sqrt(mean((mdlY - optParams.exp(idxTrial).y).^2));

            case 'simulate'
                optError=nan;
            otherwise
                assert(0,'Error: invalid optParams.name');
        end

        

        if(simConfig.flag_debugFitting==1)
    
            figure(figDebugFitting);
            subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
    
    
            switch optParams.name
                case 'responseTimeScaling'
                    txtName = expTRSS2017.activeLengtheningData(idx).seriesName;
                    i0=strfind(txtName,'Exp.');
                    txtName(1,i0:4)='Sim.';
            
                    plot(benchRecord.normFiberLength(:,idx).*lceOptMdl,...
                         benchRecord.normCrossBridgeForceAlongTendon(:,idx),...
                         '-','Color',lineColorsSimTRSS2017(idx,:),...
                         'DisplayName',...
                         txtName);
                    hold on;
                    
                case 'xeStiffnessDampingScaling'
                    txtName = expTRSS2017.activeLengtheningData(idx).seriesName;
                    i0=strfind(txtName,'Exp.');
                    txtName(1,i0:4)='Sim.';
            
                    plot(benchRecord.normFiberLength(:,idx).*lceOptMdl,...
                         benchRecord.normCrossBridgeForceAlongTendon(:,idx),...
                         '-','Color',lineColorsSimTRSS2017(idx,:),...
                         'DisplayName',...
                         txtName);
                    hold on;
                    
                case 'QToF'
                    txtName = expTRSS2017.activeLengtheningData(idx).seriesName;
                    i0=strfind(txtName,'Exp.');
                    txtName(1,i0:4)='Sim.';
            
                    plot(benchRecord.normFiberLength(:,idx).*lceOptMdl,...
                         benchRecord.normFiberForce(:,idx),...
                         '-','Color',lineColorsSimTRSS2017(idx,:),...
                         'DisplayName',...
                         txtName);
                    hold on;

                case 'QToK'
                    txtName = expTRSS2017.activeLengtheningData(idx).seriesName;
                    i0=strfind(txtName,'Exp.');
                    txtName(1,i0:4)='Sim.';
            
                    plot(benchRecord.normFiberLength(:,idx).*lceOptMdl,...
                         benchRecord.normFiberForce(:,idx),...
                         '-','Color',lineColorsSimTRSS2017(idx,:),...
                         'DisplayName',...
                         txtName);
                    hold on;
    
                    
                case 'f1HNPreload'
                    txtName = expTRSS2017.activeLengtheningData(idx).seriesName;
                    i0=strfind(txtName,'Exp.');
                    txtName(1,i0:4)='Sim.';
            
                    plot(benchRecord.normFiberLength(:,idx).*lceOptMdl,...
                         benchRecord.normFiberForce(:,idx),...
                         '-','Color',lineColorsSimTRSS2017(idx,:),...
                         'DisplayName',...
                         txtName);
                    hold on;
    
                otherwise
                    assert(0,'Error: invalid optParams.name');
            end

    
        end
    end   
