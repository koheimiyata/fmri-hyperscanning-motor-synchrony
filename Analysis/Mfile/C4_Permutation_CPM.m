%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Permutation_CPM
%
%   2026/02/06
%   Written by KMIyata
%
%   Project: InterSync
%   Purpose: Perform connectome-based predictive modeling (CPM) with
%            permutation testing to assess prediction significance.
%            Connectivity matrices from two runs are averaged, and
%            behavioral scores (circle-drawing synchrony index) are
%            predicted via leave-one-out cross-validation.
%            Edges are optionally constrained to exclude patterns
%            unique to one brain (constrain_flag).
%            Null distribution is built from 1000 permutations of
%            behavioral labels using parallel computing (parfor).
%
%   Input  : ConMatrix_fMRIPrep_Shen268_R1st/R2nd.mat,
%            SynchronyMetrics.mat
%   Output : true_prediction_r_pos/neg, pval_pos/neg, true_PM/NM,
%            Pos_rs/Neg_rs (null distributions)
%
%   Dependencies: predict_behavior.m
%
%   This script is based on 2015 Xilin Shen and Emily Finn 
%
%   Copyright 2015 Xilin Shen and Emily Finn 
%   This code is released under the terms of the GNU GPL v2. This code
%   is not FDA approved for clinical use; it is provided
%   freely for research purposes. If using this in a publication
%   please reference this properly as: 
%   Finn ES, Shen X, Scheinost D, Rosenberg MD, Huang, Chun MM,
%   Papademetris X & Constable RT. (2015). Functional connectome
%   fingerprinting: Identifying individuals using patterns of brain
%   connectivity. Nature Neuroscience 18, 1664-1671.
%   This code provides a framework for implementing functional
%   connectivity-based behavioral prediction in a leave-one-subject-out
%   cross-validation scheme, as described in Finn, Shen et al 2015 (see above
%   for full reference). The first input ('all_mats') is a pre-calculated
%   MxMxN matrix containing all individual-subject connectivity matrices,
%   where M = number of nodes in the chosen brain atlas and N = number of
%   subjects. Each element (i,j,k) in these matrices represents the
%   correlation between the BOLD timecourses of nodes i and j in subject k
%   during a single fMRI session. The second input ('all_behav') is the
%   Nx1 vector of scores for the behavior of interest for all subjects.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial settings
thresh = 0.05; 
num_roi = 268;
constrain_flag = 1; % whether the edges are all included 1:No, 0:Yes
no_iterations = 1000;
addpath(fullfile(pwd, 'Functions')) % add Functions folder to path

%% Main
direct = pwd;
% load mat file
load(fullfile(fileparts(direct), 'Mat_file', 'CPM', 'ConMatrix_fMRIPrep_Shen268_R1st.mat')) % 1st session
CM(:,:,:,1) = Zs; clear Zs

load(fullfile(fileparts(direct), 'Mat_file', 'CPM', 'ConMatrix_fMRIPrep_Shen268_R2nd.mat')) % 2nd session
CM(:,:,:,2) = Zs;

load(fullfile(fileparts(direct), 'Mat_file', 'SynchronyMetrics.mat')) % Behavioral data

% Average across two sessions
all_mats = mean(CM,4);
all_behav = mean(Synchrony.Real.SyncIndex(:,1:2), 2);

no_pairs = size(all_mats,3);

% calculate the true prediction correlation
[true_prediction_r_pos, true_prediction_r_neg, true_PM, true_NM] = predict_behavior(all_mats, all_behav, thresh, num_roi, 2, constrain_flag);

% number of iterations for permutation testing
Pos_rs = zeros(no_iterations,1);
Neg_rs = zeros(no_iterations,1);
Pos_rs(1) = true_prediction_r_pos;
Neg_rs(1) = true_prediction_r_neg;

save(fullfile('..', 'Mat_file/', 'true_PM.mat'), 'true_PM')

% create estimate distribution of the test statistic via random shuffles of data lables   
parfor it = 2:no_iterations
    fprintf('\n Performing iteration %d out of %d', it, no_iterations);
    new_behav = all_behav(randperm(no_pairs));
    [Pos_rs(it), Neg_rs(it), ~, ~] = predict_behavior(all_mats, new_behav, thresh, num_roi, 2, constrain_flag);    
end

sorted_prediction_r_pos = sort(Pos_rs,'descend');
position_pos            = find(sorted_prediction_r_pos == true_prediction_r_pos);
pval_pos                = position_pos/no_iterations;

sorted_prediction_r_neg = sort(Neg_rs,'descend');
position_neg            = find(sorted_prediction_r_neg == true_prediction_r_neg);
pval_neg                = position_neg(1)/no_iterations;
