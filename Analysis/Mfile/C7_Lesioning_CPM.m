%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C7_Lesioning_CPM
%
%   2026/04/07
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Assess the necessity of each functional network for CPM
%            prediction of interpersonal motor synchrony.
%            For each of the 10 functional networks (Shen268 parcellation),
%            all edges connected to nodes of that network — in BOTH brains
%            simultaneously — are zeroed out (computationally lesioned),
%            and CPM is re-run.
%
%            Approach mirrors Yoo et al. (2022) Fig. 5a:
%              "we computationally lesioned all the nodes in each of the
%               ten networks iteratively"
%
%            Lesioning is applied symmetrically: for functional network k,
%            nodes are removed from both subA (label k) and subB (label k+10)
%            simultaneously, preserving the inter-brain symmetry constraint.
%
%            Masking is done by zero-filling all_mats (no new matrix files).
%            Symmetry constraint (constrain_flag) is inherited from the
%            main CPM and applied automatically to remaining edges.
%
%   Performance change:
%     Δr = r_full - r_lesioned  (positive = network contributes to prediction)
%
%   Input  : ConMatrix_fMRIPrep_Shen268_R1st/R2nd.mat,
%            SynchronyMetrics.mat,
%            Shen268_N10_Labels.mat (L: 536x1, labels 1-20)
%   Output : CPM_Lesioning_results.mat
%
%   Dependencies: predict_behavior_v2.m, compute_q2.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial settings
thresh         = 0.05;
num_roi        = 268;
constrain_flag = 1;   % intersect (same as main CPM)
cmethod        = 2;   % Spearman (same as main CPM)
nNetworks      = 10;  % functional networks (subA: 1-10, subB: 11-20)
addpath(fullfile(pwd, 'Functions'))

%% Load data
direct = pwd;

load(fullfile(fileparts(direct), 'Mat_file', 'CPM', 'ConMatrix_fMRIPrep_Shen268_R1st.mat'))
CM(:,:,:,1) = Zs; clear Zs

load(fullfile(fileparts(direct), 'Mat_file', 'CPM', 'ConMatrix_fMRIPrep_Shen268_R2nd.mat'))
CM(:,:,:,2) = Zs; clear Zs

load(fullfile(fileparts(direct), 'Mat_file', 'SynchronyMetrics.mat'))

% Load network labels: L is 536x1, values 1-10 for subA, 11-20 for subB
load(fullfile(fileparts(direct), 'ROI_temp', 'Shen268_N10_Labels'))
% L(1:268) in range 1-10, L(269:536) in range 11-20

% Average across two sessions
all_mats  = mean(CM, 4); clear CM
all_behav = mean(Synchrony.Real.SyncIndex(:,1:2), 2);
no_pairs  = size(all_mats, 3);

fprintf('\n=== C7: Network Lesioning CPM ===\n');
fprintf('N pairs = %d,  N networks = %d\n', no_pairs, nNetworks);

%% Run full model (baseline)
fprintf('\n--- Full model (baseline) ---\n');
[r_pos_full, r_neg_full, ~, ~, pred_pos_full, pred_neg_full] = ...
    predict_behavior(all_mats, all_behav, thresh, num_roi, cmethod, constrain_flag);
q2_pos_full = compute_q2(pred_pos_full, all_behav);
q2_neg_full = compute_q2(pred_neg_full, all_behav);

fprintf('  Full model:  r_pos=%.4f  q2_pos=%.4f  r_neg=%.4f  q2_neg=%.4f\n', ...
        r_pos_full, q2_pos_full, r_neg_full, q2_neg_full);

%% Lesioning loop
r_pos_lesion  = zeros(nNetworks, 1);
r_neg_lesion  = zeros(nNetworks, 1);
q2_pos_lesion = zeros(nNetworks, 1);
q2_neg_lesion = zeros(nNetworks, 1);
n_lesioned    = zeros(nNetworks, 1); % number of zeroed edges per network

for k = 1:nNetworks
    fprintf('\n--- Lesioning network %d ---\n', k);

    % Identify nodes for functional network k in both brains
    subA_nodes = find(L(1:num_roi) == k);                % indices in 1:268
    subB_nodes = find(L(num_roi+1:end) == k+nNetworks) + num_roi; % indices in 269:536

    all_nodes = [subA_nodes; subB_nodes]; % combined node indices (1:536)

    % Count edges to be zeroed (each edge counted once: row or col in all_nodes)
    % Unique edges involving these nodes in the upper triangle of 536x536
    other_nodes = setdiff((1:size(all_mats,1))', all_nodes);
    n_within_lesion   = length(subA_nodes)*(length(subA_nodes)-1)/2 * 2;  % both brains
    n_between_lesion  = length(subA_nodes)^2;                              % A-k x B-k
    n_cross_lesion    = length(all_nodes) * length(other_nodes);           % k to all others
    n_lesioned(k) = n_within_lesion + n_between_lesion + n_cross_lesion;

    fprintf('  Network %d: %d subA nodes, %d subB nodes, ~%d edges zeroed\n', ...
            k, length(subA_nodes), length(subB_nodes), n_lesioned(k));

    % Zero-fill: remove all edges connected to network k nodes
    lesioned_mats = all_mats;
    lesioned_mats(all_nodes, :, :) = 0;
    lesioned_mats(:, all_nodes, :) = 0;

    % Run CPM on lesioned matrix
    [r_pos_lesion(k), r_neg_lesion(k), ~, ~, pred_pos, pred_neg] = ...
        predict_behavior_v2(lesioned_mats, all_behav, thresh, num_roi, cmethod, constrain_flag);
    q2_pos_lesion(k) = compute_q2(pred_pos, all_behav);
    q2_neg_lesion(k) = compute_q2(pred_neg, all_behav);

    fprintf('  r_pos=%.4f  q2_pos=%.4f  r_neg=%.4f  q2_neg=%.4f\n', ...
            r_pos_lesion(k), q2_pos_lesion(k), r_neg_lesion(k), q2_neg_lesion(k));
end

%% Compute performance change (Δr = full - lesioned)
delta_r_pos = r_pos_full - r_pos_lesion;
delta_r_neg = r_neg_full - r_neg_lesion;

%% Edge count control (cf. Yoo et al. Supplementary Table 3)
% Check whether Δr is driven by number of lesioned edges
[rho_pos, p_pos] = corr(n_lesioned, delta_r_pos, 'Type', 'Spearman');
[rho_neg, p_neg] = corr(n_lesioned, delta_r_neg, 'Type', 'Spearman');

fprintf('\n--- Edge count control ---\n');
fprintf('  Corr(n_lesioned_edges, Δr_pos): rho=%.3f, p=%.3f\n', rho_pos, p_pos);
fprintf('  Corr(n_lesioned_edges, Δr_neg): rho=%.3f, p=%.3f\n', rho_neg, p_neg);
if p_pos < 0.05 || p_neg < 0.05
    fprintf('  WARNING: Δr is significantly correlated with edge count.\n');
    fprintf('           Consider reporting this as a limitation.\n');
else
    fprintf('  OK: Δr is not significantly driven by edge count.\n');
end

%% Print summary table
fprintf('\n\n=== Results Summary ===\n');
fprintf('Baseline (full model):  r_pos=%.4f  q2_pos=%.4f  r_neg=%.4f  q2_neg=%.4f\n\n', ...
        r_pos_full, q2_pos_full, r_neg_full, q2_neg_full);

fprintf('%-10s  %8s  %8s  %8s  %8s  %8s  %8s\n', ...
        'Network', 'r_pos', 'q2_pos', 'r_neg', 'q2_neg', 'Δr_pos', 'Δr_neg');
fprintf('%s\n', repmat('-', 1, 72));
for k = 1:nNetworks
    fprintf('Net %2d    %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n', k, ...
            r_pos_lesion(k), q2_pos_lesion(k), r_neg_lesion(k), q2_neg_lesion(k), ...
            delta_r_pos(k), delta_r_neg(k));
end

%% Save results
results.full.r_pos          = r_pos_full;
results.full.r_neg          = r_neg_full;
results.full.q2_pos         = q2_pos_full;
results.full.q2_neg         = q2_neg_full;
results.lesion.r_pos        = r_pos_lesion;
results.lesion.r_neg        = r_neg_lesion;
results.lesion.q2_pos       = q2_pos_lesion;
results.lesion.q2_neg       = q2_neg_lesion;
results.lesion.delta_r_pos  = delta_r_pos;
results.lesion.delta_r_neg  = delta_r_neg;
results.lesion.n_lesioned   = n_lesioned;
results.edgecount_ctrl.rho_pos = rho_pos;
results.edgecount_ctrl.p_pos   = p_pos;
results.edgecount_ctrl.rho_neg = rho_neg;
results.edgecount_ctrl.p_neg   = p_neg;

out_path = fullfile(fileparts(direct), 'Mat_file', 'CPM_Lesioning_results.mat');
save(out_path, 'results', 'all_behav')
fprintf('\nResults saved: %s\n', out_path);

cd(direct)
