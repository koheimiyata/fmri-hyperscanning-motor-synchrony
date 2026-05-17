%% master.m
clear; clc;

fprintf('Starting the process...\n');

% Preprocessing
P1_fMRIPrep2SPM
P2_Discard_initialFiles
P3_SaveConfoundsAsMat

% Connectome-based predictive modeling
C1_FLevel_CPM
C2_ResidT2F
C3_ConMatrix_Shen268
C4_Permutation_CPM
C5_SummarizeMatrix

% Additional CPM analyses
C6_WithinBetween_CPM
C7_Lesioning_CPM
C8_Isolation_CPM
C9_SingleBrain_CPM

fprintf('All processes have been completed successfully.\n');
