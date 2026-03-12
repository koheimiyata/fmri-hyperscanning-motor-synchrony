function [R_pos, R_neg, PM, NM] = predict_behavior(all_mats, all_behav, thresh, num_roi, cmethod, constrain_flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   predict_behavior
%
%   2026/3/9
%   Written by KMiyata
%
%   Study: InterSync
%   Purpose: Core function for Connectome-based Predictive Modeling (CPM).
%            Implements leave-one-out cross-validation (LOOCV) to predict
%            behavioral scores from functional connectivity matrices.
%
%   Inputs:
%     all_mats       : M x M x N connectivity matrices (Fisher z-transformed)
%     all_behav      : N x 1 behavioral scores
%     thresh         : p-value threshold for edge selection (e.g., 0.05)
%     num_roi        : number of ROIs per brain (e.g., 268 for Shen atlas)
%     cmethod        : correlation method
%                      1 = robust regression, 2 = Spearman, 3 = Pearson
%     constrain_flag : edge constraint across brains
%                      0 = no constraint (all edges included)
%                      1 = intersect (keep edges significant in BOTH brains)
%                      2 = union (keep edges significant in EITHER brain)
%
%   Outputs:
%     R_pos : prediction-behavior correlation using positive edges
%     R_neg : prediction-behavior correlation using negative edges
%     PM    : M x M x N positive masks (per LOOCV fold)
%     NM    : M x M x N negative masks (per LOOCV fold)
%
%   Based on: Finn et al. (2015) Nature Neuroscience 18, 1664-1671.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2015 Xilin Shen and Emily Finn 

% This code is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes. If using this in a publication
% please reference this properly as: 

% Finn ES, Shen X, Scheinost D, Rosenberg MD, Huang, Chun MM,
% Papademetris X & Constable RT. (2015). Functional connectome
% fingerprinting: Identifying individuals using patterns of brain
% connectivity. Nature Neuroscience 18, 1664-1671.

% This code provides a framework for implementing functional
% connectivity-based behavioral prediction in a leave-one-subject-out
% cross-validation scheme, as described in Finn, Shen et al 2015 (see above
% for full reference). The first input ('all_mats') is a pre-calculated
% MxMxN matrix containing all individual-subject connectivity matrices,
% where M = number of nodes in the chosen brain atlas and N = number of
% subjects. Each element (i,j,k) in these matrices represents the
% correlation between the BOLD timecourses of nodes i and j in subject k
% during a single fMRI session. The second input ('all_behav') is the
% Nx1 vector of scores for the behavior of interest for all subjects.

% As in the reference paper, the predictive power of the model is assessed
% via correlation between predicted and observed scores across all
% subjects. Note that this assumes normal or near-normal distributions for
% both vectors, and does not assess absolute accuracy of predictions (only
% relative accuracy within the sample). It is recommended to explore
% additional/alternative metrics for assessing predictive power, such as
% prediction error sum of squares or prediction r^2.

no_sub = size(all_mats,3);
no_node = size(all_mats,1);
no_node2 = size(all_mats,2);

behav_pred_pos = zeros(no_sub,1);
behav_pred_neg = zeros(no_sub,1);
behav_pred = zeros(no_sub,1);

for leftout = 1:no_sub
    fprintf('\n Leaving out subj # %6.3f',leftout);
    
    % leave out subject from matrices and behavior
    
    train_mats = all_mats;
    train_mats(:,:,leftout) = [];
    train_vcts = reshape(train_mats,[],size(train_mats,3));
    
    train_behav = all_behav;
    train_behav(leftout) = [];
    
    if cmethod == 1
        % correlate all edges with behavior using robust regression
        edge_no = size(train_vcts,1);
        r_mat = zeros(1, edge_no);
        p_mat = zeros(1, edge_no);

        for edge_i = 1: edge_no
            [~, stats] = robustfit(train_vcts(edge_i,:)', train_behav);
            cur_t = stats.t(2);
            r_mat(edge_i) = sign(cur_t)*sqrt(cur_t^2/(no_sub-1-2+cur_t^2));
            p_mat(edge_i) = 2*(1-tcdf(abs(cur_t), no_sub-1-2));  %two tailed
        end
    
        % % correlate all edges with behavior using partial correlation
        % [r_mat, p_mat] = partialcorr(train_vcts', train_behav, age);
    
    elseif cmethod == 2
        % correlate all edges with behavior using rank correlation
        [r_mat, p_mat] = corr(train_vcts', train_behav, 'type', 'Spearman','Rows','pairwise');

    elseif cmethod == 3
        % correlate all edges with behavior using rank correlation
        [r_mat, p_mat] = corr(train_vcts', train_behav, 'Type', 'Pearson','Rows','pairwise');
    end
            
    r_mat = reshape(r_mat,no_node,[]);
    p_mat = reshape(p_mat,no_node,[]);
    
    % set threshold and define masks 
    pos_mask = zeros(no_node, no_node2);
    neg_mask = zeros(no_node, no_node2);
    
    pos_edge = find( r_mat >0 & p_mat < thresh);
    neg_edge = find( r_mat <0 & p_mat < thresh);
    
    pos_mask(pos_edge) = 1;
    neg_mask(neg_edge) = 1;
    
    if constrain_flag == 1
        % constrains due to the same condition beween two
        % within 
        for rx=1:num_roi
            for ry=1:num_roi
                if pos_mask(rx,ry)==1 && pos_mask(rx+num_roi,ry+num_roi)==0
                    pos_mask(rx,ry)=0;
                elseif pos_mask(rx,ry)==0 && pos_mask(rx+num_roi,ry+num_roi)==1
                    pos_mask(rx+num_roi,ry+num_roi)=0;
                end
                if neg_mask(rx,ry)==1 && neg_mask(rx+num_roi,ry+num_roi)==0
                   neg_mask(rx,ry)=0;
                elseif neg_mask(rx,ry)==0 && neg_mask(rx+num_roi,ry+num_roi)==1
                   neg_mask(rx+num_roi,ry+num_roi)=0;
                end
            end
        end
        % between
        for rx=1:num_roi
            for ry=num_roi+1:num_roi*2
                if pos_mask(rx,ry)==1 && pos_mask(rx+num_roi,ry-num_roi)==0
                    pos_mask(rx,ry)=0;
                    pos_mask(ry,rx)=0;
                elseif pos_mask(rx,ry)==0 && pos_mask(rx+num_roi,ry-num_roi)==1               
                    pos_mask(rx+num_roi,ry-num_roi)=0;
                    pos_mask(ry-num_roi,rx+num_roi)=0;
                end
                if neg_mask(rx,ry)==1 && neg_mask(rx+num_roi,ry-num_roi)==0
                    neg_mask(rx,ry)=0;
                    neg_mask(ry,rx)=0; 
                elseif neg_mask(rx,ry)==0 && neg_mask(rx+num_roi,ry-num_roi)==1                
                    neg_mask(rx+num_roi,ry-num_roi)=0;
                    neg_mask(ry-num_roi,rx+num_roi)=0;
                end
            end
        end
    elseif constrain_flag == 2
        % constrains due to the same condition beween two
        % within 
        for rx=1:num_roi
            for ry=1:num_roi
                if pos_mask(rx,ry)==1 && pos_mask(rx+num_roi,ry+num_roi)==0
                    pos_mask(rx+num_roi,ry+num_roi)=1;
                elseif pos_mask(rx,ry)==0 && pos_mask(rx+num_roi,ry+num_roi)==1
                    pos_mask(rx,ry)=1;
                end
                if neg_mask(rx,ry)==1 && neg_mask(rx+num_roi,ry+num_roi)==0
                   neg_mask(rx+num_roi,ry+num_roi)=1;
                elseif neg_mask(rx,ry)==0 && neg_mask(rx+num_roi,ry+num_roi)==1
                   neg_mask(rx,ry)=1;
                end
            end
        end
        % between
        for rx=1:num_roi
            for ry=num_roi+1:num_roi*2
                if pos_mask(rx,ry)==1 && pos_mask(rx+num_roi,ry-num_roi)==0
                    pos_mask(rx+num_roi,ry-num_roi)=1;
                    pos_mask(ry-num_roi,rx+num_roi)=1;
                elseif pos_mask(rx,ry)==0 && pos_mask(rx+num_roi,ry-num_roi)==1               
                    pos_mask(rx,ry)=1;
                    pos_mask(ry,rx)=1;
                end
                if neg_mask(rx,ry)==1 && neg_mask(rx+num_roi,ry-num_roi)==0
                    neg_mask(rx+num_roi,ry-num_roi)=1;
                    neg_mask(ry-num_roi,rx+num_roi)=1;
                elseif neg_mask(rx,ry)==0 && neg_mask(rx+num_roi,ry-num_roi)==1                
                    neg_mask(rx,ry)=1;
                    neg_mask(ry,rx)=1;
                end
            end
        end

    end

    % get sum of all edges in TRAIN subs (divide by 2 to control for the
    % fact that matrices are symmetric)
    train_sumpos = zeros(no_sub-1,1);
    train_sumneg = zeros(no_sub-1,1);
    
    for ss = 1:size(train_sumpos,1)
        train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask,"omitnan"),"omitnan")/2;
        train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask,"omitnan"),"omitnan")/2;
    end
    
    % run model on TEST sub   
    test_mat = all_mats(:,:,leftout);
     
    % behav_pred(leftout) = b(1)*test_sumpos + b(2)*test_sumneg + b(3);
    if sum(sum(pos_mask))~=0
        fit_pos = polyfit(train_sumpos, train_behav, 1); % build model on TRAIN subs
        test_sumpos = sum(sum(test_mat.*pos_mask,"omitnan"),"omitnan")/2;
        behav_pred_pos(leftout) = fit_pos(1)*test_sumpos + fit_pos(2);
    else
        behav_pred_pos(leftout) = NaN;
    end
    if sum(sum(neg_mask))~=0
        fit_neg = polyfit(train_sumneg, train_behav, 1); % build model on TRAIN subs
        test_sumneg = sum(sum(test_mat.*neg_mask,"omitnan"),"omitnan")/2;
        behav_pred_neg(leftout) = fit_neg(1)*test_sumneg + fit_neg(2);
    else
        behav_pred_neg(leftout) = NaN;
    end

    PM(:,:,leftout)=pos_mask;
    NM(:,:,leftout)=neg_mask;
end

% compare predicted and observed scores
if sum(isnan(behav_pred_pos))<=3
    [R_pos, ~] = corr(behav_pred_pos,all_behav,'rows','complete');
else
    R_pos=0;
end
if sum(isnan(behav_pred_neg))<=3
    [R_neg, ~] = corr(behav_pred_neg,all_behav,'rows','complete');
else
    R_neg=0;
end