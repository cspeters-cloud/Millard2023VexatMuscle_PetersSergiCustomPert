%% PS_runSingleAgonistPert_HumanSoleus

%%
disp("Running Single Agonist Perturbations with Human Soleus")

%% wrapper inputs
muscleArchitecture = defaultHumanSoleus.musculotendon ;
sarcomereProperties = defaultHumanSoleus.sarcomere ;
normMuscleCurves = defaultHumanSoleus.curves ;

%%
act = PS_inputState.activationState.a ;
dact = PS_inputState.activationState.da ;

lp = [PS_inputState.pathState.lp_stretch PS_inputState.pathState.lp_shorten] ;
dlp = [PS_inputState.pathState.dlp_stretch PS_inputState.pathState.dlp_shorten] ;

%% initialize state
modelConfig0 = struct( ...
              'iterMax'                 , 100             , ...
              'tol'                     , 1e-3   , ... 
              'tolInit'                 , sqrt(eps)       , ...
              'minActivation'           , 0.0             , ...
              'useElasticTendon'        , flag_useElasticTendon , ...
              'initializeState'         , 1                );          
              activationState0 = [dact(1); act(1)] ;
              pathState0 = [dlp(1, 1); lp(1, 1)] ;
              if flag_useElasticTendon
                  muscleState0 = [0;0;0;0] ;
              else
                  muscleState0 = [0;0;0] ;
              end
              
mtInfo0 = calcMillard2023VexatMuscleInfo(activationState0,...
                                                      pathState0,...
                                                      muscleState0,...
                                                      muscleArchitecture,...
                                                      sarcomereProperties,...
                                                      normMuscleCurves,...
                                                      modelConfig0);
                                                  
%% state integrator loop
modelConfig1 = struct( ...
              'iterMax'                 , 100             , ...
              'tol'                     , 1e-3   , ... 
              'tolInit'                 , sqrt(eps)       , ...
              'minActivation'           , 0.0             , ...
              'useElasticTendon'        , flag_useElasticTendon , ...
              'initializeState'         , 0               );          
              
mtInfo1 = cell(1001, 2);
    
for ip = 1 : 2

    mtInfo1{1, ip} = mtInfo0;

    activationState1 = [dact(2); act(2)] ;
    pathState1 = [dlp(2,ip); lp(2,ip)] ;
    if flag_useElasticTendon
        muscleState1 = [mtInfo0.state.value(1); mtInfo0.state.value(2); mtInfo0.state.value(3); mtInfo0.state.value(4)] ;
    else
        muscleState1 = [mtInfo0.state.value(1); mtInfo0.state.value(2); mtInfo0.state.value(3)] ;
    end

    for it = 2 : 1001

        mtInfo1{it, ip} = calcMillard2023VexatMuscleInfo(activationState1,...
                                                          pathState1,...
                                                          muscleState1,...
                                                          muscleArchitecture,...
                                                          sarcomereProperties,...
                                                          normMuscleCurves,...
                                                          modelConfig1);
        activationState1 = [dact(it); act(it)];
        pathState1 = [dlp(it, ip); lp(it, ip)] ;
        
        if flag_useElasticTendon
            muscleState1 = [mtInfo1{it, ip}.state.value(1); ...
            mtInfo1{it, ip}.state.value(2); ...
            mtInfo1{it, ip}.state.value(3); ...
            mtInfo1{it, ip}.state.value(4)] ;
        else
            muscleState1 = [mtInfo1{it, ip}.state.value(1); ...
            mtInfo1{it, ip}.state.value(2); ...
            mtInfo1{it, ip}.state.value(3)] ;
        end

    end    
    
end

%% reshape to timeseries
mtInfoTimeSeries.muscleLengthInfo = struct() ;
    mtInfoTimeSeries.muscleLengthInfo.tendonLength = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.normTendonLength = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.tendonStrain = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.normFiberSlidingLength = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.fiberSlidingLength = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.fiberLength = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.fiberLengthAlongTendon = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.normFiberLength = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.pennationAngle = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.cosPennationAngle = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.sinPennationAngle = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.titin1Length = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.normTitin1Length = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.titin2Length = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.normTitin2Length = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.crossBridgeLength = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.normCrossBridgeLength = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.isClamped = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.fiberPassiveForceLengthMultiplier = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.fiberActiveForceLengthMultiplier = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.normFiberSlidingVelocity = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.fiberSlidingVelocity = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.crossBridgeVelocity = zeros(1001,2);
    mtInfoTimeSeries.muscleLengthInfo.normCrossBridgeVelocity = zeros(1001,2);

mtInfoTimeSeries.fiberVelocityInfo = struct() ;
    mtInfoTimeSeries.fiberVelocityInfo.tendonVelocity = zeros(1001,2);
    mtInfoTimeSeries.fiberVelocityInfo.normTendonVelocity = zeros(1001,2);
    mtInfoTimeSeries.fiberVelocityInfo.fiberVelocity = zeros(1001,2);
    mtInfoTimeSeries.fiberVelocityInfo.fiberVelocityAlongTendon = zeros(1001,2);
    mtInfoTimeSeries.fiberVelocityInfo.normFiberVelocity = zeros(1001,2);
    mtInfoTimeSeries.fiberVelocityInfo.pennationAngularVelocity = zeros(1001,2);
    mtInfoTimeSeries.fiberVelocityInfo.fiberForceVelocityMultiplier = zeros(1001,2);
    mtInfoTimeSeries.fiberVelocityInfo.titin1Velocity = zeros(1001,2);
    mtInfoTimeSeries.fiberVelocityInfo.normTitin1Velocity = zeros(1001,2);
    mtInfoTimeSeries.fiberVelocityInfo.titin2Velocity = zeros(1001,2);
    mtInfoTimeSeries.fiberVelocityInfo.normTitin2Velocity = zeros(1001,2);

mtInfoTimeSeries.muscleDynamicsInfo = struct() ;
    mtInfoTimeSeries.muscleDynamicsInfo.activation = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.activationDerivative = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.lambda = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.dlambda = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.tendonForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normTendonForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.tendonStiffness = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.tendonDamping = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.fiberForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.fiberForceAlongTendon = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normFiberForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.activeFiberForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normActiveFiberForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.passiveFiberForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normPassiveFiberForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.crossBridgeStiffness = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.crossBridgeDamping = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.fiberActivePower = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.fiberPassivePower = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.tendonPower = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.musclePower = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.fiberParallelElementPower = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.ecmStiffness = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.ecmDamping = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.titin1Stiffness = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.titin1Damping = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.titin2Stiffness = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.titin2Damping = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.titin1Force = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normTitin1Force = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.titin2Force = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normTitin2Force = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.titin1SpringForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normTitin1SpringForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.titin2SpringForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normTitin2SpringForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.titin1DampingForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normTitin1DampingForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.titin2DampingForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normTitin2DampingForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.crossBridgeForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normCrossBridgeForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.crossBridgeSpringForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normCrossBridgeSpringForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.crossBridgeDampingForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normCrossBridgeDampingForce = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.dampingForces = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.dampingPower = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.boundaryPower = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.fiberStiffness = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.fiberDamping = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normFiberStiffness = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normFiberDamping = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.musculotendonStiffness = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.musculotendonDamping = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normMusculotendonStiffness = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.normMusculotendonDamping = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.dT = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.dW = zeros(1001,2);
    mtInfoTimeSeries.muscleDynamicsInfo.dV = zeros(1001,2);

mtInfoTimeSeries.state = struct() ;
    if flag_useElasticTendon
        mtInfoTimeSeries.state.value = zeros(1001,2,4);
        mtInfoTimeSeries.state.derivative = zeros(1001,2,4);
    else
        mtInfoTimeSeries.state.value = zeros(1001,2,3);
        mtInfoTimeSeries.state.derivative = zeros(1001,2,3);
    end
    
for it = 1 : 1001
    
    for ip = 1 : 2
        
        mtInfoTimeSeries.muscleLengthInfo.tendonLength(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.tendonLength;
        mtInfoTimeSeries.muscleLengthInfo.normTendonLength(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.normTendonLength;
        mtInfoTimeSeries.muscleLengthInfo.tendonStrain(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.tendonStrain;
        mtInfoTimeSeries.muscleLengthInfo.normFiberSlidingLength(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.normFiberSlidingLength;
        mtInfoTimeSeries.muscleLengthInfo.fiberSlidingLength(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.fiberSlidingLength;
        mtInfoTimeSeries.muscleLengthInfo.fiberLength(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.fiberLength;
        mtInfoTimeSeries.muscleLengthInfo.fiberLengthAlongTendon(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.fiberLengthAlongTendon;
        mtInfoTimeSeries.muscleLengthInfo.normFiberLength(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.normFiberLength;
        mtInfoTimeSeries.muscleLengthInfo.pennationAngle(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.pennationAngle;
        mtInfoTimeSeries.muscleLengthInfo.cosPennationAngle(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.cosPennationAngle;
        mtInfoTimeSeries.muscleLengthInfo.sinPennationAngle(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.sinPennationAngle;
        mtInfoTimeSeries.muscleLengthInfo.titin1Length(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.titin1Length;
        mtInfoTimeSeries.muscleLengthInfo.normTitin1Length(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.normTitin1Length;
        mtInfoTimeSeries.muscleLengthInfo.titin2Length(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.titin2Length;
        mtInfoTimeSeries.muscleLengthInfo.normTitin2Length(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.normTitin2Length;
        mtInfoTimeSeries.muscleLengthInfo.crossBridgeLength(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.crossBridgeLength;
        mtInfoTimeSeries.muscleLengthInfo.normCrossBridgeLength(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.normCrossBridgeLength;
        mtInfoTimeSeries.muscleLengthInfo.isClamped(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.isClamped;
        mtInfoTimeSeries.muscleLengthInfo.fiberPassiveForceLengthMultiplier(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.fiberPassiveForceLengthMultiplier;
        mtInfoTimeSeries.muscleLengthInfo.fiberActiveForceLengthMultiplier(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.fiberActiveForceLengthMultiplier;
        mtInfoTimeSeries.muscleLengthInfo.normFiberSlidingVelocity(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.normFiberSlidingVelocity;
        mtInfoTimeSeries.muscleLengthInfo.fiberSlidingVelocity(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.fiberSlidingVelocity;
        mtInfoTimeSeries.muscleLengthInfo.crossBridgeVelocity(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.crossBridgeVelocity;
        mtInfoTimeSeries.muscleLengthInfo.normCrossBridgeVelocity(it,ip) = mtInfo1{it,ip}.muscleLengthInfo.normCrossBridgeVelocity;

        mtInfoTimeSeries.fiberVelocityInfo.tendonVelocity(it,ip) = mtInfo1{it,ip}.fiberVelocityInfo.tendonVelocity;
        mtInfoTimeSeries.fiberVelocityInfo.normTendonVelocity(it,ip) = mtInfo1{it,ip}.fiberVelocityInfo.normTendonVelocity;
        mtInfoTimeSeries.fiberVelocityInfo.fiberVelocity(it,ip) = mtInfo1{it,ip}.fiberVelocityInfo.fiberVelocity;
        mtInfoTimeSeries.fiberVelocityInfo.fiberVelocityAlongTendon(it,ip) = mtInfo1{it,ip}.fiberVelocityInfo.fiberVelocityAlongTendon;
        mtInfoTimeSeries.fiberVelocityInfo.normFiberVelocity(it,ip) = mtInfo1{it,ip}.fiberVelocityInfo.normFiberVelocity;
        mtInfoTimeSeries.fiberVelocityInfo.pennationAngularVelocity(it,ip) = mtInfo1{it,ip}.fiberVelocityInfo.pennationAngularVelocity;
        mtInfoTimeSeries.fiberVelocityInfo.fiberForceVelocityMultiplier(it,ip) = mtInfo1{it,ip}.fiberVelocityInfo.fiberForceVelocityMultiplier;
        mtInfoTimeSeries.fiberVelocityInfo.titin1Velocity(it,ip) = mtInfo1{it,ip}.fiberVelocityInfo.titin1Velocity;
        mtInfoTimeSeries.fiberVelocityInfo.normTitin1Velocity(it,ip) = mtInfo1{it,ip}.fiberVelocityInfo.normTitin1Velocity;
        mtInfoTimeSeries.fiberVelocityInfo.titin2Velocity(it,ip) = mtInfo1{it,ip}.fiberVelocityInfo.titin2Velocity;
        mtInfoTimeSeries.fiberVelocityInfo.normTitin2Velocity(it,ip) = mtInfo1{it,ip}.fiberVelocityInfo.normTitin2Velocity;

        mtInfoTimeSeries.muscleDynamicsInfo.activation(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.activation;
        mtInfoTimeSeries.muscleDynamicsInfo.activationDerivative(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.activationDerivative;
        mtInfoTimeSeries.muscleDynamicsInfo.lambda(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.lambda;
        mtInfoTimeSeries.muscleDynamicsInfo.dlambda(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.dlambda;
        mtInfoTimeSeries.muscleDynamicsInfo.tendonForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.tendonForce;
        mtInfoTimeSeries.muscleDynamicsInfo.normTendonForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normTendonForce;
        mtInfoTimeSeries.muscleDynamicsInfo.tendonStiffness(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.tendonStiffness;
        mtInfoTimeSeries.muscleDynamicsInfo.tendonDamping(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.tendonDamping;
        mtInfoTimeSeries.muscleDynamicsInfo.fiberForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.fiberForce;
        mtInfoTimeSeries.muscleDynamicsInfo.fiberForceAlongTendon(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.fiberForceAlongTendon;
        mtInfoTimeSeries.muscleDynamicsInfo.normFiberForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normFiberForce;
        mtInfoTimeSeries.muscleDynamicsInfo.activeFiberForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.activeFiberForce;
        mtInfoTimeSeries.muscleDynamicsInfo.normActiveFiberForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normActiveFiberForce;
        mtInfoTimeSeries.muscleDynamicsInfo.passiveFiberForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.passiveFiberForce;
        mtInfoTimeSeries.muscleDynamicsInfo.normPassiveFiberForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normPassiveFiberForce;
        mtInfoTimeSeries.muscleDynamicsInfo.crossBridgeStiffness(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.crossBridgeStiffness;
        mtInfoTimeSeries.muscleDynamicsInfo.crossBridgeDamping(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.crossBridgeDamping;
        mtInfoTimeSeries.muscleDynamicsInfo.fiberActivePower(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.fiberActivePower;
        mtInfoTimeSeries.muscleDynamicsInfo.fiberPassivePower(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.fiberPassivePower;
        mtInfoTimeSeries.muscleDynamicsInfo.tendonPower(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.tendonPower;
        mtInfoTimeSeries.muscleDynamicsInfo.musclePower(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.musclePower;
        mtInfoTimeSeries.muscleDynamicsInfo.fiberParallelElementPower(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.fiberParallelElementPower;
        mtInfoTimeSeries.muscleDynamicsInfo.ecmStiffness(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.ecmStiffness;
        mtInfoTimeSeries.muscleDynamicsInfo.ecmDamping(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.ecmDamping;
        mtInfoTimeSeries.muscleDynamicsInfo.titin1Stiffness(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.titin1Stiffness;
        mtInfoTimeSeries.muscleDynamicsInfo.titin1Damping(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.titin1Damping;
        mtInfoTimeSeries.muscleDynamicsInfo.titin2Stiffness(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.titin2Stiffness;
        mtInfoTimeSeries.muscleDynamicsInfo.titin2Damping(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.titin2Damping;
        mtInfoTimeSeries.muscleDynamicsInfo.titin1Force(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.titin1Force;
        mtInfoTimeSeries.muscleDynamicsInfo.normTitin1Force(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normTitin1Force;
        mtInfoTimeSeries.muscleDynamicsInfo.titin2Force(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.titin2Force;
        mtInfoTimeSeries.muscleDynamicsInfo.normTitin2Force(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normTitin2Force;
        mtInfoTimeSeries.muscleDynamicsInfo.titin1SpringForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.titin1SpringForce;
        mtInfoTimeSeries.muscleDynamicsInfo.normTitin1SpringForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normTitin1SpringForce;
        mtInfoTimeSeries.muscleDynamicsInfo.titin2SpringForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.titin2SpringForce;
        mtInfoTimeSeries.muscleDynamicsInfo.normTitin2SpringForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normTitin2SpringForce;
        mtInfoTimeSeries.muscleDynamicsInfo.titin1DampingForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.titin1DampingForce;
        mtInfoTimeSeries.muscleDynamicsInfo.normTitin1DampingForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normTitin1DampingForce;
        mtInfoTimeSeries.muscleDynamicsInfo.titin2DampingForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.titin2DampingForce;
        mtInfoTimeSeries.muscleDynamicsInfo.normTitin2DampingForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normTitin2DampingForce;
        mtInfoTimeSeries.muscleDynamicsInfo.crossBridgeForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.crossBridgeForce;
        mtInfoTimeSeries.muscleDynamicsInfo.normCrossBridgeForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normCrossBridgeForce;
        mtInfoTimeSeries.muscleDynamicsInfo.crossBridgeSpringForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.crossBridgeSpringForce;
        mtInfoTimeSeries.muscleDynamicsInfo.normCrossBridgeSpringForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normCrossBridgeSpringForce;
        mtInfoTimeSeries.muscleDynamicsInfo.crossBridgeDampingForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.crossBridgeDampingForce;
        mtInfoTimeSeries.muscleDynamicsInfo.normCrossBridgeDampingForce(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normCrossBridgeDampingForce;
        mtInfoTimeSeries.muscleDynamicsInfo.dampingForces(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.dampingForces;
        mtInfoTimeSeries.muscleDynamicsInfo.dampingPower(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.dampingPower;
        mtInfoTimeSeries.muscleDynamicsInfo.boundaryPower(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.boundaryPower;
        mtInfoTimeSeries.muscleDynamicsInfo.fiberStiffness(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.fiberStiffness;
        mtInfoTimeSeries.muscleDynamicsInfo.fiberDamping(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.fiberDamping;
        mtInfoTimeSeries.muscleDynamicsInfo.normFiberStiffness(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normFiberStiffness;
        mtInfoTimeSeries.muscleDynamicsInfo.normFiberDamping(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normFiberDamping;
        mtInfoTimeSeries.muscleDynamicsInfo.musculotendonStiffness(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.musculotendonStiffness;
        mtInfoTimeSeries.muscleDynamicsInfo.musculotendonDamping(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.musculotendonDamping;
        mtInfoTimeSeries.muscleDynamicsInfo.normMusculotendonStiffness(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normMusculotendonStiffness;
        mtInfoTimeSeries.muscleDynamicsInfo.normMusculotendonDamping(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.normMusculotendonDamping;
        mtInfoTimeSeries.muscleDynamicsInfo.dT(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.dT;
        mtInfoTimeSeries.muscleDynamicsInfo.dW(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.dW;
        mtInfoTimeSeries.muscleDynamicsInfo.dV(it,ip) = mtInfo1{it,ip}.muscleDynamicsInfo.dV;
        
        mtInfoTimeSeries.state.value(it,ip,:) = mtInfo1{it,ip}.state.value;
        mtInfoTimeSeries.state.derivative(it,ip,:) = mtInfo1{it,ip}.state.derivative;

        
    end
end

