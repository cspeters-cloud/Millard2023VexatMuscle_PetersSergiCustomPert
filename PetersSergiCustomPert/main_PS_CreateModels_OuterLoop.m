%% main_PS_CreateModels_OuterLoop_forMillard

%%
% rootDir         = getRootProjectDirectory();
% projectFolders  = getProjectFolders(rootDir);

%% global model parameters
flag_useOctave            = 0; 

flag_makeAndSavePubPlots  = 0;
% plotOutputFolder          = [projectFolders.output_plots_MuscleCurves,filesep];

normMaxActiveTitinToActinDamping = 65;

normPevkToActinAttachmentPointDefault = 0.5;

normFiberLengthAtOneNormPassiveForceDefault = 1.367732948060934e+00;

normFiberLengthAtOneNormPassiveForceDefault_RT = 1.421393259171109e+00;

ecmForceFractionDefault = 0.56;

ecmForceFractionRabbitPsoas     = 1-(0.728-0.158); %Fig 8A, Prado et al.
ecmForceFractionHumanSoleus     = ecmForceFractionDefault;
ecmForceFractionFelineSoleus    = ecmForceFractionDefault;

smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));
flag_enableNumericallyNonZeroGradients  = 1;

wlcTitinModel = 1;

linearTitinModel=0;

useCalibratedCurves     = 1;

useTwoSidedTitinCurves  = 0;

%% human soleus model parameters
flag_plotAllHumanSoleusCurves           = 0;
scaleOptimalFiberLengthHumanSoleus      = 1.0; 
scaleMaximumIsometricTensionHumanSoleus = 1;

%% paths
addpath( genpath(projectFolders.parameters)     );
addpath( genpath(projectFolders.curves)         );
addpath( genpath(projectFolders.experiments)    );
addpath( genpath(projectFolders.simulation)     );
addpath( genpath(projectFolders.models)         );
addpath( genpath(projectFolders.postprocessing) );

%% plot configuration
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

%% human soleus and achilles tendon model
fprintf('\n\nCreating: default human soleus model\n');
fprintf('\tused to simulate the Ig and PEVK kinematics from Trombitas et al.\n\n');

defaultHumanSoleus = createHumanSoleusModel(...
                        normPevkToActinAttachmentPointDefault,...
                        normMaxActiveTitinToActinDamping,...                        
                        normFiberLengthAtOneNormPassiveForceDefault,... 
                        ecmForceFractionHumanSoleus,...
                        linearTitinModel,...
                        useCalibratedCurves,...
                        useTwoSidedTitinCurves,...
                        smallNumericallyNonZeroNumber,...
                        flag_enableNumericallyNonZeroGradients,...
                        scaleOptimalFiberLengthHumanSoleus,...
                        scaleMaximumIsometricTensionHumanSoleus,...
                        projectFolders,...
                        flag_useOctave);
                    
defaultHumanSoleus.musculotendon.tendonSlackLength = 0.1 * defaultHumanSoleus.musculotendon.tendonSlackLength ;

filePathHumanSoleus = fullfile(projectFolders.output_structs_FittedModels,...
                                'defaultHumanSoleus.mat');
save(filePathHumanSoleus,'defaultHumanSoleus');

if(flag_plotAllHumanSoleusCurves==1)
    figHumanSoleusCurves = ...
    plotStructOfBezierSplines( defaultHumanSoleus.curves,...
                                      {'Inverse','use'});       
end

%% plotting
if(flag_makeAndSavePubPlots==1)
  [success] = plotMuscleCurves( fittedFelineSoleusKBR1994Fig12_ET,...
                                defaultHumanSoleus,...
                                activeForceLengthCurveAnnotationPoints,...
                                felineSoleusActiveForceLengthDataDefault,...
                                felineSoleusPassiveForceLengthDataDefault,...
                                normFiberLengthAtOneNormPassiveForceDefault,...
                                plotOutputFolder,...
                                projectFolders);
end 
