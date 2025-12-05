%% main_PS_OuterLoop

%%
clc ;
close all ;
clear all ;

%% paths
rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);
addpath( genpath(projectFolders.parameters)     ) ;
addpath( genpath(projectFolders.curves)         ) ;
addpath( genpath(projectFolders.experiments)    ) ;
addpath( genpath(projectFolders.simulation)     ) ;
addpath( genpath(projectFolders.models)         ) ;
addpath( genpath(projectFolders.postprocessing) ) ;

addpath([rootDir '\PetersSergiCustomPert']) ;

%% create models
main_PS_CreateModels_OuterLoop ;

%% single agonist perturbations: human soleus
% select_useElasticTendon = 1 ; % 1 for Elastic Tendon
select_useElasticTendon = 0 ; % 0 for Rigid Tendon
main_PS_SingleAgonistPert_HumanSoleus_OuterLoop ;

%%
close all ;
PS_plotSingleHSolPerturbation ;
