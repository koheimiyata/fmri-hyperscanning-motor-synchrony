function q2 = compute_q2(behav_pred, all_behav)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   compute_q2
%
%   2026/04/07
%   Written by KMiyata
%
%   Purpose: Compute prediction q² (cross-validated R²)
%
%            q² = 1 - SS_residual / SS_total
%
%            where SS_residual = Σ(y_pred - y_obs)²
%                  SS_total    = Σ(y_obs  - mean(y_obs))²
%
%            Values < 0 (model worse than predicting the mean) are set to 0.
%
%   Inputs:
%     behav_pred : N x 1 predicted behavioral scores (may contain NaN)
%     all_behav  : N x 1 observed behavioral scores
%
%   Output:
%     q2 : scalar prediction q²  (0 if model is below-chance)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

valid = ~isnan(behav_pred);

% Need at least 3 valid predictions
if sum(valid) < 3
    q2 = 0;
    return;
end

y_pred = behav_pred(valid);
y_obs  = all_behav(valid);

ss_res = sum((y_pred - y_obs).^2);
ss_tot = sum((y_obs  - mean(y_obs)).^2);

if ss_tot == 0
    q2 = 0;
    return;
end

q2 = 1 - ss_res / ss_tot;

% Clip negative values to 0 (below-chance performance)
if q2 < 0
    q2 = 0;
end
