function [R_pos, R_neg, PM, NM, behav_pred_pos, behav_pred_neg] = predict_behavior(all_mats, all_behav, thresh, num_roi, cmethod, constrain_flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   predict_behavior
%
%   2026/04/07
%   Modified by KMiyata
%
%   Study: InterSync
%   Purpose: Extended version of predict_behavior.m.
%            Returns predicted behavioral scores (behav_pred_pos,
%            behav_pred_neg) in addition to the original outputs,
%            enabling computation of prediction q² (cross-validated R²)
%            as used in Yoo et al. (2022) Nature Human Behaviour.
%
%   Inputs:  (same as predict_behavior.m)
%     all_mats       : M x M x N connectivity matrices (Fisher z-transformed)
%     all_behav      : N x 1 behavioral scores
%     thresh         : p-value threshold for edge selection (e.g., 0.05)
%     num_roi        : number of ROIs per brain (e.g., 268 for Shen atlas)
%     cmethod        : 1 = robust regression, 2 = Spearman, 3 = Pearson
%     constrain_flag : 0 = no constraint, 1 = intersect, 2 = union
%
%   Outputs:
%     R_pos          : prediction-behavior correlation (positive edges)
%     R_neg          : prediction-behavior correlation (negative edges)
%     PM             : M x M x N positive masks (per LOOCV fold)
%     NM             : M x M x N negative masks (per LOOCV fold)
%     behav_pred_pos : N x 1 predicted scores (positive network)
%     behav_pred_neg : N x 1 predicted scores (negative network)
%
%   Note: Use compute_q2.m to compute q² from behav_pred and all_behav.
%
%   Based on: Finn et al. (2015) Nature Neuroscience 18, 1664-1671.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

no_sub  = size(all_mats, 3);
no_node = size(all_mats, 1);
no_node2 = size(all_mats, 2);

behav_pred_pos = zeros(no_sub, 1);
behav_pred_neg = zeros(no_sub, 1);

for leftout = 1:no_sub
    fprintf('\n Leaving out subj # %6.3f', leftout);

    % Leave out one pair from matrices and behavior
    train_mats = all_mats;
    train_mats(:,:,leftout) = [];
    train_vcts = reshape(train_mats, [], size(train_mats, 3));

    train_behav = all_behav;
    train_behav(leftout) = [];

    % Correlate edges with behavior
    if cmethod == 1
        edge_no = size(train_vcts, 1);
        r_mat = zeros(1, edge_no);
        p_mat = zeros(1, edge_no);
        for edge_i = 1:edge_no
            [~, stats] = robustfit(train_vcts(edge_i,:)', train_behav);
            cur_t = stats.t(2);
            r_mat(edge_i) = sign(cur_t) * sqrt(cur_t^2 / (no_sub-1-2 + cur_t^2));
            p_mat(edge_i) = 2 * (1 - tcdf(abs(cur_t), no_sub-1-2));
        end
    elseif cmethod == 2
        [r_mat, p_mat] = corr(train_vcts', train_behav, 'type', 'Spearman', 'Rows', 'pairwise');
    elseif cmethod == 3
        [r_mat, p_mat] = corr(train_vcts', train_behav, 'Type', 'Pearson', 'Rows', 'pairwise');
    end

    r_mat = reshape(r_mat, no_node, []);
    p_mat = reshape(p_mat, no_node, []);

    % Define edge masks
    pos_mask = zeros(no_node, no_node2);
    neg_mask = zeros(no_node, no_node2);

    pos_edge = find(r_mat > 0 & p_mat < thresh);
    neg_edge = find(r_mat < 0 & p_mat < thresh);

    pos_mask(pos_edge) = 1;
    neg_mask(neg_edge) = 1;

    % Apply symmetry constraint
    if constrain_flag == 1
        % Intersect: keep edges significant in BOTH brains
        % Within-brain
        for rx = 1:num_roi
            for ry = 1:num_roi
                if pos_mask(rx,ry)==1 && pos_mask(rx+num_roi,ry+num_roi)==0
                    pos_mask(rx,ry) = 0;
                elseif pos_mask(rx,ry)==0 && pos_mask(rx+num_roi,ry+num_roi)==1
                    pos_mask(rx+num_roi,ry+num_roi) = 0;
                end
                if neg_mask(rx,ry)==1 && neg_mask(rx+num_roi,ry+num_roi)==0
                    neg_mask(rx,ry) = 0;
                elseif neg_mask(rx,ry)==0 && neg_mask(rx+num_roi,ry+num_roi)==1
                    neg_mask(rx+num_roi,ry+num_roi) = 0;
                end
            end
        end
        % Between-brain
        for rx = 1:num_roi
            for ry = num_roi+1:num_roi*2
                if pos_mask(rx,ry)==1 && pos_mask(rx+num_roi,ry-num_roi)==0
                    pos_mask(rx,ry) = 0;
                    pos_mask(ry,rx) = 0;
                elseif pos_mask(rx,ry)==0 && pos_mask(rx+num_roi,ry-num_roi)==1
                    pos_mask(rx+num_roi,ry-num_roi) = 0;
                    pos_mask(ry-num_roi,rx+num_roi) = 0;
                end
                if neg_mask(rx,ry)==1 && neg_mask(rx+num_roi,ry-num_roi)==0
                    neg_mask(rx,ry) = 0;
                    neg_mask(ry,rx) = 0;
                elseif neg_mask(rx,ry)==0 && neg_mask(rx+num_roi,ry-num_roi)==1
                    neg_mask(rx+num_roi,ry-num_roi) = 0;
                    neg_mask(ry-num_roi,rx+num_roi) = 0;
                end
            end
        end

    elseif constrain_flag == 2
        % Union: keep edges significant in EITHER brain
        % Within-brain
        for rx = 1:num_roi
            for ry = 1:num_roi
                if pos_mask(rx,ry)==1 && pos_mask(rx+num_roi,ry+num_roi)==0
                    pos_mask(rx+num_roi,ry+num_roi) = 1;
                elseif pos_mask(rx,ry)==0 && pos_mask(rx+num_roi,ry+num_roi)==1
                    pos_mask(rx,ry) = 1;
                end
                if neg_mask(rx,ry)==1 && neg_mask(rx+num_roi,ry+num_roi)==0
                    neg_mask(rx+num_roi,ry+num_roi) = 1;
                elseif neg_mask(rx,ry)==0 && neg_mask(rx+num_roi,ry+num_roi)==1
                    neg_mask(rx,ry) = 1;
                end
            end
        end
        % Between-brain
        for rx = 1:num_roi
            for ry = num_roi+1:num_roi*2
                if pos_mask(rx,ry)==1 && pos_mask(rx+num_roi,ry-num_roi)==0
                    pos_mask(rx+num_roi,ry-num_roi) = 1;
                    pos_mask(ry-num_roi,rx+num_roi) = 1;
                elseif pos_mask(rx,ry)==0 && pos_mask(rx+num_roi,ry-num_roi)==1
                    pos_mask(rx,ry) = 1;
                    pos_mask(ry,rx) = 1;
                end
                if neg_mask(rx,ry)==1 && neg_mask(rx+num_roi,ry-num_roi)==0
                    neg_mask(rx+num_roi,ry-num_roi) = 1;
                    neg_mask(ry-num_roi,rx+num_roi) = 1;
                elseif neg_mask(rx,ry)==0 && neg_mask(rx+num_roi,ry-num_roi)==1
                    neg_mask(rx,ry) = 1;
                    neg_mask(ry,rx) = 1;
                end
            end
        end
    end

    % Compute edge strength sums for training pairs
    train_sumpos = zeros(no_sub-1, 1);
    train_sumneg = zeros(no_sub-1, 1);
    for ss = 1:size(train_sumpos, 1)
        train_sumpos(ss) = sum(sum(train_mats(:,:,ss) .* pos_mask, 'omitnan'), 'omitnan') / 2;
        train_sumneg(ss) = sum(sum(train_mats(:,:,ss) .* neg_mask, 'omitnan'), 'omitnan') / 2;
    end

    % Predict test pair
    test_mat = all_mats(:,:,leftout);

    if sum(sum(pos_mask)) ~= 0
        fit_pos = polyfit(train_sumpos, train_behav, 1);
        test_sumpos = sum(sum(test_mat .* pos_mask, 'omitnan'), 'omitnan') / 2;
        behav_pred_pos(leftout) = fit_pos(1) * test_sumpos + fit_pos(2);
    else
        behav_pred_pos(leftout) = NaN;
    end
    if sum(sum(neg_mask)) ~= 0
        fit_neg = polyfit(train_sumneg, train_behav, 1);
        test_sumneg = sum(sum(test_mat .* neg_mask, 'omitnan'), 'omitnan') / 2;
        behav_pred_neg(leftout) = fit_neg(1) * test_sumneg + fit_neg(2);
    else
        behav_pred_neg(leftout) = NaN;
    end

    PM(:,:,leftout) = pos_mask;
    NM(:,:,leftout) = neg_mask;
end

% Prediction-behavior correlation
if sum(isnan(behav_pred_pos)) <= 3
    [R_pos, ~] = corr(behav_pred_pos, all_behav, 'rows', 'complete');
else
    R_pos = 0;
end
if sum(isnan(behav_pred_neg)) <= 3
    [R_neg, ~] = corr(behav_pred_neg, all_behav, 'rows', 'complete');
else
    R_neg = 0;
end
