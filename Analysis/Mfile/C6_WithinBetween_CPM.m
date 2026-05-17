%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C6_WithinBetween_CPM
%
%   2026/04/07
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Test whether a two-brain network is necessary for predicting
%            interpersonal motor synchrony by comparing CPM performance
%            using (1) the full connectivity matrix, (2) within-brain
%            edges only, and (3) between-brain edges only.
%
%            Approach: The original 536x536 connectivity matrix is masked
%            by zero-filling irrelevant blocks — no new connectivity
%            matrices are created (memory efficient).
%
%            Matrix block structure (536 = 268 subA + 268 subB):
%              [  within-A  |  between  ]    rows/cols  1-268   = subA
%              [  between'  |  within-B ]    rows/cols 269-536  = subB
%
%            Conditions:
%              Full         : all edges (baseline)
%              Within-only  : zero out off-diagonal (between) blocks
%              Between-only : zero out diagonal (within) blocks
%
%            Note on edge counts (before symmetry constraint):
%              Within-brain : 2 x 268x267/2 = 71,556 edges
%              Between-brain:   268x268      = 71,824 edges
%              Difference is only 268 (~0.4%), so edge count is not a
%              meaningful confound for this particular comparison.
%
%   Performance metrics:
%     - Pearson's r  : correlation between predicted and observed scores
%     - Prediction q²: cross-validated R² (negative values set to 0)
%
%   Note: Permutation testing is not conducted because the small sample
%         size in LOPO-CV limits the resolution of the null distribution.
%         Results are reported descriptively.
%
%   Input  : ConMatrix_fMRIPrep_Shen268_R1st/R2nd.mat,
%            SynchronyMetrics.mat
%   Output : CPM_WithinBetween_results.mat
%            true_PM_full.mat, true_PM_within.mat, true_PM_between.mat
%            (PM files are compatible with C5_SummarizeMatrix.m)
%
%   Dependencies: predict_behavior_v2.m, compute_q2.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial settings
thresh         = 0.05;
num_roi        = 268;
constrain_flag = 1;   % 1: intersect constraint (same as main CPM)
cmethod        = 2;   % 2: Spearman (same as main CPM)
addpath(fullfile(pwd, 'Functions'))

%% Load data
direct = pwd;

load(fullfile(fileparts(direct), 'Mat_file', 'CPM', 'ConMatrix_fMRIPrep_Shen268_R1st.mat'))
CM(:,:,:,1) = Zs; clear Zs

load(fullfile(fileparts(direct), 'Mat_file', 'CPM', 'ConMatrix_fMRIPrep_Shen268_R2nd.mat'))
CM(:,:,:,2) = Zs; clear Zs

load(fullfile(fileparts(direct), 'Mat_file', 'SynchronyMetrics.mat'))

% Average across two sessions
all_mats   = mean(CM, 4); clear CM
all_behav  = mean(Synchrony.Real.SyncIndex(:,1:2), 2);
no_pairs   = size(all_mats, 3);
no_node    = size(all_mats, 1); % 536

fprintf('\n=== C6: Within / Between Brain CPM ===\n');
fprintf('N pairs = %d\n', no_pairs);

%% Create masked matrices (zero-fill — no new matrix files created)

% (1) Full matrix: no masking
all_mats_full = all_mats;

% (2) Within-brain only: zero out between-brain blocks
%     Preserve within-SubA [1:268, 1:268] and within-SubB [269:536, 269:536]
all_mats_within = all_mats;
all_mats_within(1:num_roi, num_roi+1:end, :) = 0;  % top-right (A→B)
all_mats_within(num_roi+1:end, 1:num_roi, :) = 0;  % bottom-left (B→A)

% (3) Between-brain only: zero out within-brain blocks
%     Preserve between-brain blocks [1:268, 269:536] and [269:536, 1:268]
all_mats_between = all_mats;
all_mats_between(1:num_roi, 1:num_roi, :)               = 0;  % within-SubA
all_mats_between(num_roi+1:end, num_roi+1:end, :)       = 0;  % within-SubB

%% Edge count report (before symmetry constraint & thresholding)
n = num_roi;
n_within  = 2 * n*(n-1)/2;   % both brains
n_between = n * n;            % full A×B block (not upper-tri; matrix is asymmetric within block)
fprintf('\nEdge counts (before constraint):\n');
fprintf('  Within-brain : %d\n', n_within);
fprintf('  Between-brain: %d\n', n_between);
fprintf('  Difference   : %d (%.1f%%)\n', abs(n_within - n_between), ...
        100*abs(n_within - n_between)/n_between);

%% Run CPM for each condition
conditions  = {'Full', 'Within-only', 'Between-only'};
mats_list   = {all_mats_full, all_mats_within, all_mats_between};
results     = struct();

for ci = 1:3
    fprintf('\n--- Condition: %s ---\n', conditions{ci});

    [r_pos, r_neg, PM, NM, pred_pos, pred_neg] = ...
        predict_behavior(mats_list{ci}, all_behav, thresh, num_roi, cmethod, constrain_flag);

    q2_pos = compute_q2(pred_pos, all_behav);
    q2_neg = compute_q2(pred_neg, all_behav);

    % Store results
    cname = strrep(conditions{ci}, '-', '_');
    cname = strrep(cname, ' ', '_');
    results.(cname).r_pos      = r_pos;
    results.(cname).r_neg      = r_neg;
    results.(cname).q2_pos     = q2_pos;
    results.(cname).q2_neg     = q2_neg;
    results.(cname).pred_pos   = pred_pos;
    results.(cname).pred_neg   = pred_neg;

    % Save PM (compatible with C5_SummarizeMatrix.m)
    true_PM = PM;
    pm_name = fullfile(fileparts(direct), 'Mat_file', ...
                       ['true_PM_' lower(strrep(conditions{ci},' ','_')) '.mat']);
    save(pm_name, 'true_PM')
    fprintf('  PM saved: %s\n', pm_name);
end

%% Print summary table
fprintf('\n\n=== Results Summary ===\n');
fprintf('%-18s  %8s  %8s  %8s  %8s\n', 'Condition', 'r_pos', 'q2_pos', 'r_neg', 'q2_neg');
fprintf('%s\n', repmat('-', 1, 60));
for ci = 1:3
    cname = strrep(conditions{ci}, '-', '_');
    cname = strrep(cname, ' ', '_');
    fprintf('%-18s  %8.4f  %8.4f  %8.4f  %8.4f\n', conditions{ci}, ...
            results.(cname).r_pos, results.(cname).q2_pos, ...
            results.(cname).r_neg, results.(cname).q2_neg);
end

%% Save results
out_path = fullfile(fileparts(direct), 'Mat_file', 'CPM_WithinBetween_results.mat');
save(out_path, 'results', 'conditions', 'all_behav')
fprintf('\nResults saved: %s\n', out_path);
fprintf('\nNote: PM files saved as true_PM_full/within/between.mat.\n');
fprintf('      Run C5_SummarizeMatrix.m after changing the load line to inspect\n');
fprintf('      network-level edge contributions for each condition.\n');

cd(direct)
