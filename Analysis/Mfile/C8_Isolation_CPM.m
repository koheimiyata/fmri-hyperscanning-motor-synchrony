%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C8_Isolation_CPM
%
%   2026/04/07
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Assess the sufficiency of each functional network for CPM
%            prediction of interpersonal motor synchrony.
%
%            For each of 10 functional networks, CPM is run using ONLY
%            edges where both endpoints belong to that network
%            (within-network subset; cf. Yoo et al. 2022, Fig. 5b).
%
%            For functional network k, the subset includes:
%              - within-SubA-k  : edges among SubA nodes of network k
%              - within-SubB-k  : edges among SubB nodes of network k
%              - between-brain-k: edges between SubA_k and SubB_k
%                                 (unique to hyperscanning)
%
%            Masking is done by zero-filling all_mats (no new matrix files).
%            Symmetry constraint (constrain_flag) is inherited from main CPM.
%
%   Edge count note:
%     n_edges(k) = nk*(2*nk - 1)  where nk = nodes per brain in network k
%       = nk*(nk-1)/2 [within-SubA]
%       + nk*(nk-1)/2 [within-SubB]
%       + nk^2        [between-brain, k-to-k]
%
%   Performance metrics: Pearson's r and prediction q².
%   Note: Permutation testing not conducted (see C4 for rationale).
%         Results are reported descriptively.
%
%   Input  : ConMatrix_fMRIPrep_Shen268_R1st/R2nd.mat,
%            SynchronyMetrics.mat,
%            Shen268_N10_Labels.mat (L: 536x1, labels 1-20)
%   Output : CPM_Subset_results.mat
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
nNetworks      = 10;
addpath(fullfile(pwd, 'Functions'))

%% Load data
direct = pwd;

load(fullfile(fileparts(direct), 'Mat_file', 'CPM', 'ConMatrix_fMRIPrep_Shen268_R1st.mat'))
CM(:,:,:,1) = Zs; clear Zs

load(fullfile(fileparts(direct), 'Mat_file', 'CPM', 'ConMatrix_fMRIPrep_Shen268_R2nd.mat'))
CM(:,:,:,2) = Zs; clear Zs

load(fullfile(fileparts(direct), 'Mat_file', 'SynchronyMetrics.mat'))
load(fullfile(fileparts(direct), 'ROI_temp', 'Shen268_N10_Labels'))
% L: 536x1, values 1-10 for SubA nodes, 11-20 for SubB nodes

all_mats  = mean(CM, 4); clear CM
all_behav = mean(Synchrony.Real.SyncIndex(:,1:2), 2);
no_pairs  = size(all_mats, 3);
no_node   = size(all_mats, 1); % 536

fprintf('\n=== C8: Network Subset CPM (Type A: within-network) ===\n');
fprintf('N pairs = %d,  N networks = %d\n', no_pairs, nNetworks);

%% Preallocate
r_pos   = zeros(nNetworks, 1);
r_neg   = zeros(nNetworks, 1);
q2_pos  = zeros(nNetworks, 1);
q2_neg  = zeros(nNetworks, 1);
n_edges = zeros(nNetworks, 1);

%% Subset loop
for k = 1:nNetworks
    fprintf('\n--- Network %d ---\n', k);

    % Node indices for functional network k
    subA_k = find(L(1:num_roi) == k);                          % in 1:268
    subB_k = find(L(num_roi+1:end) == k+nNetworks) + num_roi;  % in 269:536
    nk = length(subA_k);

    % Edge count (undirected):
    %   within-SubA + within-SubB + between-brain A_k x B_k
    n_edges(k) = 2 * nk*(nk-1)/2 + nk^2;  % = nk*(2*nk - 1)

    % Zero-fill: keep only edges within network k (both endpoints in k)
    subset_mats = zeros(no_node, no_node, no_pairs);
    subset_mats(subA_k, subA_k, :) = all_mats(subA_k, subA_k, :);  % within-SubA-k
    subset_mats(subB_k, subB_k, :) = all_mats(subB_k, subB_k, :);  % within-SubB-k
    subset_mats(subA_k, subB_k, :) = all_mats(subA_k, subB_k, :);  % between-brain k→k
    subset_mats(subB_k, subA_k, :) = all_mats(subB_k, subA_k, :);  % symmetric

    % Run CPM
    [r_pos(k), r_neg(k), ~, ~, pred_pos, pred_neg] = ...
        predict_behavior(subset_mats, all_behav, thresh, num_roi, cmethod, constrain_flag);
    q2_pos(k) = compute_q2(pred_pos, all_behav);
    q2_neg(k) = compute_q2(pred_neg, all_behav);

    fprintf('  r_pos=%6.4f  q2_pos=%6.4f  r_neg=%6.4f  q2_neg=%6.4f  (n_edges=%d, nk=%d)\n', ...
            r_pos(k), q2_pos(k), r_neg(k), q2_neg(k), n_edges(k), nk);
end

%% Edge count control
[rho_pos, p_pos] = corr(n_edges, r_pos, 'Type', 'Spearman');
[rho_neg, p_neg] = corr(n_edges, r_neg, 'Type', 'Spearman');

fprintf('\n--- Edge count control ---\n');
fprintf('  Corr(n_edges, r_pos): rho=%.3f p=%.3f\n', rho_pos, p_pos);
fprintf('  Corr(n_edges, r_neg): rho=%.3f p=%.3f\n', rho_neg, p_neg);
if p_pos < 0.05 || p_neg < 0.05
    fprintf('  WARNING: r correlates with edge count. Report in limitations.\n');
else
    fprintf('  OK: r is not significantly driven by edge count.\n');
end

%% Print summary table
fprintf('\n\n=== Results Summary ===\n');
fprintf('%-10s  %8s  %8s  %8s  %8s  %8s\n', ...
        'Network', 'r_pos', 'q2_pos', 'r_neg', 'q2_neg', 'n_edges');
fprintf('%s\n', repmat('-', 1, 64));
for k = 1:nNetworks
    fprintf('Net %2d    %8.4f  %8.4f  %8.4f  %8.4f  %8d\n', k, ...
            r_pos(k), q2_pos(k), r_neg(k), q2_neg(k), n_edges(k));
end

%% Save results
results.r_pos   = r_pos;
results.r_neg   = r_neg;
results.q2_pos  = q2_pos;
results.q2_neg  = q2_neg;
results.n_edges = n_edges;
results.edgecount_ctrl = struct('rho_pos', rho_pos, 'p_pos', p_pos, ...
                                'rho_neg', rho_neg, 'p_neg', p_neg);

out_path = fullfile(fileparts(direct), 'Mat_file', 'CPM_Subset_results.mat');
save(out_path, 'results', 'all_behav')
fprintf('\nResults saved: %s\n', out_path);

cd(direct)
