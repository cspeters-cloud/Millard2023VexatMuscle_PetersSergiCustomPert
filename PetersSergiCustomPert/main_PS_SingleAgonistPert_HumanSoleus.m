%% main_PS_SingleAgonistPert_HumanSoleus

%%
flag_outerLoopMode = 1;

% rootDir         = getRootProjectDirectory();
% projectFolders  = getProjectFolders(rootDir);


if(flag_outerLoopMode == 0)
  clc;
  close all;
  clear all;
  rootDir         = getRootProjectDirectory();
  projectFolders  = getProjectFolders(rootDir);

  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                      = 1;
  
  flag_fitToFig3KirchBoskovRymer1994            = 0; 
  flag_useElasticTendon                         = select_useElasticTendon;
  flag_useFiberDamping                          = 1;
  
  flag_useSameRandomPerturbationAsPublication   = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;
  flag_useCalibratedVexatCurves                = 1;
  flag_useTwoSidedTitinCurves                   = 1;

  flag_generateRandomInput    = 0;
  flag_processInputFunctions  = 0;
end 

%%
modelStructsDir = [projectFolders.output_structs_FittedModels,filesep];

addpath( genpath(projectFolders.parameters)     );
addpath( genpath(projectFolders.curves)         );
addpath( genpath(projectFolders.experiments)    );
addpath( genpath(projectFolders.simulation)     );
addpath( genpath(projectFolders.models)         );
addpath( genpath(projectFolders.postprocessing) );

%% generate perturbation waveforms
% simulation setup
sim_time = 0 : 0.001 : 1 ;
pert_mag_m = 0.0016 ; % see Kirsch, Boskov, Rymer 1994
response_dur = 0.04 ; % see Kirsch, Boskov, Rymer 1994
response_dur_ms = response_dur * 1000 ;
pert_time = sim_time(1 : response_dur_ms);
onset = 0.5 ;
onset_ms = onset * 1000;
l_M_0 = defaultHumanSoleus.musculotendon.optimalFiberLength ;
l_P_0 = l_M_0 * cos(defaultHumanSoleus.musculotendon.pennationAngle) + defaultHumanSoleus.musculotendon.tendonSlackLength ;

% inputs
load("ps_inputState.mat") ;

PS_runSingleAgonistPert_HumanSoleus;