%% main_PS_SingleAgonistPert_HumanSoleus_OuterLoop

%%
flag_OuterOuterLoopMode =1;
if(flag_OuterOuterLoopMode ==0)
    clc;
    close all;
    clear all;
end

disp('----------------------------------------');
disp(' running main_PS_SingleAgonistPert_HumanSoleus_OuterLoop_forMillard');
disp('----------------------------------------');

% rootDir         = getRootProjectDirectory();
% projectFolders  = getProjectFolders(rootDir);

%% Run all of the simulations

flag_runSimulations     = 1;

flag_useCalibratedVexatCurves = 1;
flag_useTwoSidedTitinCurves   = 0;

if(flag_runSimulations == 1) 

  %Proposed model only
  %Rigid tendon  
  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                      = 1;
  flag_fitToFig3KirchBoskovRymer1994            = 0;
  flag_useElasticTendon                         = select_useElasticTendon;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  flag_generateRandomInput    = 0;
  flag_processInputFunctions  = 0;

  main_PS_SingleAgonistPert_HumanSoleus;

  close all;
  
end