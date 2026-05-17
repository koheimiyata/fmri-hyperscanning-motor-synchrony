%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C9_SingleBrain_CPM
%
%   2026/04/07
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Assess the added value of hyperscanning (two-brain network)
%            by running CPM using only a single brain's connectivity matrix
%            per pair.
%
%            For each pair, either SubA [rows/cols 1:268] or SubB
%            [rows/cols 269:536] is selected at random (independently
%            per pair). CPM is then run on the resulting 268x268 x N
%            single-brain connectivity matrices.
%
%            To remove dependence on a single arbitrary A/B selection,
%            this procedure is repeated no_iter times (default: 1000).
%            Mean and SD of r and q² across iterations are reported.
%
%            Comparison with the two-brain CPM (C4/C6) quantifies how
%            much predictive power is gained by including inter-brain
%            connectivity.
%
%   Key differences from two-brain CPM:
%     - Matrix size: 268x268 (not 536x536)
%     - No symmetry constraint (constrain_flag = 0)
%     - Single-brain edge selection only
%
%   Performance metrics: Pearson's r and prediction q².
%
%   Input  : ConMatrix_fMRIPrep_Shen268_R1st/R2nd.mat,
%            SynchronyMetrics.mat
%   Output : CPM_SingleBrain_results.mat
%
%   Dependencies: predict_behavior_v2.m, compute_q2.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial settings
thresh         = 0.05;
num_roi        = 268;
constrain_flag = 0;    % no symmetry constraint (single brain)
cmethod        = 2;    % Spearman (same as main CPM)
no_iter        = 1000; % number of random A/B selection iterations
rng('shuffle')         % random seed (for reproducibility, replace with rng(N))
addpath(fullfile(pwd, 'Functions'))

%% Load data
direct = pwd;

load(fullfile(fileparts(direct), 'Mat_file', 'CPM', 'ConMatrix_fMRIPrep_Shen268_R1st.mat'))
CM(:,:,:,1) = Zs; clear Zs

load(fullfile(fileparts(direct), 'Mat_file', 'CPM', 'ConMatrix_fMRIPrep_Shen268_R2nd.mat'))
CM(:,:,:,2) = Zs; clear Zs

load(fullfile(fileparts(direct), 'Mat_file', 'SynchronyMetrics.mat'))

all_mats  = mean(CM, 4); clear CM
all_behav = mean(Synchrony.Real.SyncIndex(:,1:2), 2);
no_pairs  = size(all_mats, 3);

fprintf('\n=== C9: Single-Brain CPM (N=%d iterations) ===\n', no_iter);
fprintf('N pairs = %d\n', no_pairs);

%% Separate SubA and SubB matrices in advance
%  SubA: rows/cols 1:268, SubB: rows/cols 269:536
mats_A = all_mats(1:num_roi,          1:num_roi,          :); % 268x268xN
mats_B = all_mats(num_roi+1:num_roi*2, num_roi+1:num_roi*2, :); % 268x268xN

%% Preallocate result arrays
r_pos_iter  = zeros(no_iter, 1);
r_neg_iter  = zeros(no_iter, 1);
q2_pos_iter = zeros(no_iter, 1);
q2_neg_iter = zeros(no_iter, 1);

%% Main loop (parfor for speed)
fprintf('Running %d iterations...\n', no_iter);

parfor it = 1:no_iter
    % Randomly select SubA or SubB for each pair (independently)
    use_B = rand(no_pairs, 1) > 0.5; % logical: true = use SubB

    single_brain_mats = zeros(num_roi, num_roi, no_pairs);
    for ps = 1:no_pairs
        if use_B(ps)
            single_brain_mats(:,:,ps) = mats_B(:,:,ps);
        else
            single_brain_mats(:,:,ps) = mats_A(:,:,ps);
        end
    end

    % Run CPM (no symmetry constraint)
    [r_p, r_n, ~, ~, pred_pos, pred_neg] = predict_behavior(single_brain_mats, all_behav, thresh, num_roi, cmethod, constrain_flag);

    r_pos_iter(it)  = r_p;
    r_neg_iter(it)  = r_n;
    q2_pos_iter(it) = compute_q2(pred_pos, all_behav);
    q2_neg_iter(it) = compute_q2(pred_neg, all_behav);
end

%% Summarize across iterations
mean_r_pos  = mean(r_pos_iter);
mean_r_neg  = mean(r_neg_iter);
mean_q2_pos = mean(q2_pos_iter);
mean_q2_neg = mean(q2_neg_iter);
sd_r_pos    = std(r_pos_iter);
sd_r_neg    = std(r_neg_iter);
sd_q2_pos   = std(q2_pos_iter);
sd_q2_neg   = std(q2_neg_iter);

%% Print results
fprintf('\n=== Results (mean ± SD across %d iterations) ===\n', no_iter);
fprintf('%-12s  %16s  %16s\n', '', 'Positive network', 'Negative network');
fprintf('%-12s  %8s %7s  %8s %7s\n', '', 'r', 'q²', 'r', 'q²');
fprintf('%s\n', repmat('-', 1, 56));
fprintf('%-12s  %8.4f %7.4f  %8.4f %7.4f\n', ...
        'Mean', mean_r_pos, mean_q2_pos, mean_r_neg, mean_q2_neg);
fprintf('%-12s  %8.4f %7.4f  %8.4f %7.4f\n', ...
        'SD',   sd_r_pos,   sd_q2_pos,   sd_r_neg,   sd_q2_neg);

%% Save results
results.r_pos_iter  = r_pos_iter;
results.r_neg_iter  = r_neg_iter;
results.q2_pos_iter = q2_pos_iter;
results.q2_neg_iter = q2_neg_iter;
results.mean_r_pos  = mean_r_pos;
results.mean_r_neg  = mean_r_neg;
results.mean_q2_pos = mean_q2_pos;
results.mean_q2_neg = mean_q2_neg;
results.sd_r_pos    = sd_r_pos;
results.sd_r_neg    = sd_r_neg;
results.sd_q2_pos   = sd_q2_pos;
results.sd_q2_neg   = sd_q2_neg;
results.no_iter     = no_iter;

out_path = fullfile(fileparts(direct), 'Mat_file', 'CPM_SingleBrain_results.mat');
save(out_path, 'results', 'all_behav')
fprintf('\nResults saved: %s\n', out_path);

cd(direct)
