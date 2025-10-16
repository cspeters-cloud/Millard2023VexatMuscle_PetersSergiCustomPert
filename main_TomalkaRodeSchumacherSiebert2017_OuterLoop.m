clc; 
close all;
clear all;

%
% Uses the fast time constant for shortening and lengthening.
% There is no transient, but the XE follows fl*fv
%6

simConfigInput.runFitting              = 0; 
simConfigInput.generatePlots           = 1;
simConfigInput.fitToIndividualTrials   = 1; 
simConfigInput.manuallySetTimeConstant = 0;

main_TomalkaRodeSchumacherSiebert2017;

clc;
close all;
clear all;

simConfigInput.runFitting              = 1; 
simConfigInput.generatePlots           = 1;
simConfigInput.fitToIndividualTrials   = 0; 
simConfigInput.manuallySetTimeConstant = 0;

main_TomalkaRodeSchumacherSiebert2017;

%
% Uses separate time constants for shortening and lengthening.
% There is a transient, but the XE does not follows fl*fv
%

% simConfigInput.runFitting              = 1; 
% simConfigInput.generatePlots           = 1;
% simConfigInput.fitToIndividualTrials   = 1; 
% simConfigInput.manuallySetTimeConstant = 1;
% 
% main_TomalkaRodeSchumacherSiebert2017;
% 
% clc;
% close all;
% clear all;
% 
% simConfigInput.runFitting              = 1; 
% simConfigInput.generatePlots           = 1;
% simConfigInput.fitToIndividualTrials   = 0; 
% simConfigInput.manuallySetTimeConstant = 1;
% 
% main_TomalkaRodeSchumacherSiebert2017;
