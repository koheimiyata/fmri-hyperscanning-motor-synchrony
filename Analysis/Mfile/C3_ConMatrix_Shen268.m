%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C3_ConMatrix_Shen268
%
%   2026/02/06
%   Written by KMIyata
%
%   Project: InterSync
%   Purpose: Compute intra- and inter-brain ROI-to-ROI connectivity matrices
%            for each pair using the Shen 268-ROI parcellation.
%            For each ROI, the time series is summarized by first eigenvariate
%            (PCA) or mean across voxels (controlled by PCA_mean flag).
%            SubA and SubB time series are concatenated (536 ROIs total) and
%            a full correlation matrix is computed, Fisher z-transformed,
%            and saved per run.
%
%   Input  : 4D NIfTI files in Group/InterBrain/CPM/<run>/,
%            Shen268 mask (rShen_268_parcellation.nii)
%   Output : ConMatrix_fMRIPrep_Shen268_<run>.mat (Zs: 536x536xN)
%
%   Dependencies: first_eigenvariate.m, r_to_zfisher.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial settings
runs = {'R1st', 'R2nd'};
PCA_mean = 1; % 1: PCA, 2: Mean
num_regions = 268; % Number of regions

%% Main
direct = pwd;
% Load Shen268 map
mask = niftiread(fullfile('..', 'ROI_temp', 'rShen_268_parcellation.nii'));

% Mask dimension is reduced from 4D to 2D
% mask2: (X*Y*Z) x 1 vector; voxel-wise ROI label
mask2 = reshape(mask, size(mask,1)*size(mask,2)*size(mask,3), size(mask,4));

mkdir(fullfile('..', 'Mat_file', 'CPM'))

for rn = 1:2 % runs
    Zs = [];
    cd(fullfile(fileparts(direct), 'Group', 'InterBrain', 'CPM', runs{rn})) 
   
    % Data is divided into VA and VB
    va = dir("P*_subA_*.nii");
    vb = dir("P*_subB_*.nii");

    % Calcurate number of subjects
    num_pairs = size(va,1);
        
    for ps = 1:num_pairs % pairs
        % Read VA's data
        VA = va(ps).name; % e.g. P02_subA_R1st_R4D.nii
        pair_id = VA(1:3); % e.g. 'P02'
        a = double(niftiread(VA)); % Read 4D data of nth individual
        a2 = reshape(a, size(a,1)*size(a,2)*size(a,3), size(a,4)); % Reshape from 4D to 2D
        
        % Find the matching VB file explicitly and read it
        vb_match = find(strcmp({vb.name}, [pair_id '_subB_' runs{rn} '_R4D.nii']));
        VB = vb(vb_match).name;
        b = double(niftiread(VB)); % Read 4D data of nth individual
        b2 = reshape(b, size(b,1)*size(b,2)*size(b,3), size(b,4)); % Reshape from 4D to 2D
        
        % Make variable for mean or first eigenvariate of each ROI
        a2_roi_summary = zeros(size(a,4), num_regions);
        b2_roi_summary = zeros(size(b,4), num_regions);
        
        fprintf("ps=%03d  VA=%s  VB=%s\n", ps, VA, VB);
        % Calculate correlation in each ROI
        for roi_num = 1:num_regions % rois
            
            % Pick up one ROI
            mask_roi = (mask2 == roi_num);
            
            % a2 and b2 store picked up ROI
            roia2 = a2(mask_roi == 1, :);
            roib2 = b2(mask_roi == 1, :);
            
            % Summarize data (mean or PCA) default is PCA
            if PCA_mean == 1
                roia2_summary = first_eigenvariate(roia2');
                roib2_summary = first_eigenvariate(roib2');
            elseif PCA_mean == 2
                roia2_summary = mean(roia2, 1, 'omitnan');
                roib2_summary = mean(roib2, 1, 'omitnan');
            end
            
            % roia2_summary, roib2_summary 367 (time points) x 268 (ROIs) double
            if PCA_mean == 1
                a2_roi_summary(:, roi_num) = roia2_summary;
                b2_roi_summary(:, roi_num) = roib2_summary;
            elseif PCA_mean == 2
                a2_roi_summary(:, roi_num) = roia2_summary';
                b2_roi_summary(:, roi_num) = roib2_summary';
            end      
        end % rois
        
        % Conjunction a2_roi_summary and b2_roi_summary
        % roi_summary is t x 536 (268 + 268) time points x ROI*2
        roi_summary_data = [a2_roi_summary b2_roi_summary];
        
        % Calculate all correlations in a matrix
        % r_values, p_values are 536 double
        [r_values, ~] = corr(roi_summary_data); % Calculate Peason's correlation
        z_values = r_to_zfisher(r_values); % Fisher's r to z transform
        Zs(:,:,ps) = z_values; % storing z values into Zs

        clear z_values roi_mean r_values a a2 a2_roi_summary b b2 b2_roi_summary
    end
    
    matName = fullfile(fileparts(direct), 'Mat_file', 'CPM', ['ConMatrix_fMRIPrep_Shen268_' runs{rn} '.mat']);
    save(matName, 'Zs')
end % runs
cd(direct)