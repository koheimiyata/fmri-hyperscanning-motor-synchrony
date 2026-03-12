%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SummarizeMatrix
%
%   2026/3/9 
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Summarize CPM positive-mask edges at the canonical network
%            level using the Shen 268 parcellation (20 networks).
%            For each LOOCV fold, selected edges are scaled by the ratio
%            of model fraction to structural fraction (S), then averaged
%            across folds (meanS). Edge stability across folds is also
%            computed; network-pair contributions are thresholded by
%            stability (0.75, 0.80, 0.90, 1.00) and saved as CSV files.
%
%   NOTE:
%   - S(i,j,r) = [C(i,j,r)/totalE(r)] / [T(i,j)/totalPairs]
%     (fold-wise scaling, then mean over folds)
%   - Significance is based on stability across folds (S >= 1).
%
%   Input  : true_PM.mat (536x536x27 binary),
%            Shen268_N10_Labels.mat (L: 536x1, labels 1-20)
%   Output : wMatrix.csv, bMatrix.csv,
%            ContributionMatrix_<thresh>.csv (thresh = 1.00)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial settings
n = 536;      % number of nodes
nL = 20;      % number of network labels
nIter = 27;   % number of LOOCV folds
thrsh = 1;

%% Main
% Load data
load('../Mat_file/true_PM')            % PM: 536x536x27 binary matrices
load('../ROI_temp/Shen268_N10_Labels') % L: 536x1, integer label (1 to 20)
PM = true_PM;

totalPairs = n*(n-1)/2;  % total possible undirected edges
%% Preallocate
C = zeros(nL, nL, nIter);  % edge counts per network-pair, per fold
T = zeros(nL, nL);         % total possible edges per network-pair (observed labels)
totalE = zeros(nIter,1);   % total number of selected edges per fold

%% Main loop (observed)
for i = 1:n
    for j = i+1:n  % avoid double counting
        li = L(i);
        lj = L(j);

        % structural total under observed labels
        T(li, lj) = T(li, lj) + 1;
        T(lj, li) = T(lj, li) + 1;

        for r = 1:nIter
            if PM(i, j, r) == 1
                C(li, lj, r) = C(li, lj, r) + 1;
                C(lj, li, r) = C(lj, li, r) + 1;
                totalE(r) = totalE(r) + 1;
            end
        end
    end
end

%% Scaling: compute S per fold, then mean over folds (same as original)
S = zeros(nL, nL, nIter);
for r = 1:nIter
    for i = 1:nL
        for j = 1:nL
            if T(i,j) > 0 && totalE(r) > 0
                model_frac = C(i,j,r) / totalE(r);
                struct_frac = T(i,j) / totalPairs;
                S(i,j,r) = model_frac / struct_frac;
            else
                S(i,j,r) = 0;
            end
        end
    end
end

meanS = mean(S, 3);
meanC = mean(C, 3);

%% Visualization
figure(1)
heatmap(tril(meanS), 'Colormap', hot);
title('Scaled Edge Ratio (S)')

within  = meanS(1:10,1:10);
between = meanS(11:20,1:10);

writematrix(within,'wMatrix.csv')
writematrix(between,'bMatrix.csv')

c_stability = sum(S >= 1, 3);
stability = c_stability/nIter;

sig_S = stability >= thrsh;
sig_contri = meanS.*sig_S;
writematrix(sig_contri,['ContributionMatrix_' num2str(thrsh) '.csv'])