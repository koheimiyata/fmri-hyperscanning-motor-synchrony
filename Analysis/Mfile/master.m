%% master.m
clear; clc;

fprintf('Starting the process...\n');

% Preprocessing
P1_fMRIPrep2SPM
P2_Discard_initialFiles
P3_SaveConfoundsAsMat
P4_Smooth

% Behavioral regressors
B1_RelativePhase
B2_HRFConvolution

% Connectome-based predictive modeling
C1_FLevel_CPM
C2_ResidT2F
C3_ConMatrix_Shen268
C4_Permutation_CPM
C5_SummarizeMatrix

% GLM analysis
G1_FLevel_GLM
G2_CManager
G3_SLevel_GLM

fprintf('All processes have been completed successfully.\n');