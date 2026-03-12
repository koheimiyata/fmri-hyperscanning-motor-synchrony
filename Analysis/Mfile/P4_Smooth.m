%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   P4_Smooth
%
%   2026/02/06
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Apply spatial smoothing to preprocessed fMRI volumes using SPM.
%            All 3D NIfTI files (filenames starting with 'sub') across both
%            runs are collected per subject and smoothed in a single batch.
%            A Gaussian kernel with the specified FWHM (mm) is applied.
%            Smoothed files are prefixed with 's' + FWHM size (e.g., 's5').
%
%   Input  : 3D NIfTI volumes in each run directory
%   Output : Smoothed NIfTI volumes (e.g., s5*.nii), SPM batch .mat file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 

%% Initial setting
FWHMsize = 5;
p_end = 27; % the number of pairs
pairs = {'P02','P03','P04','P05','P06','P07','P08','P09', ...
    'P10','P11','P12','P13','P14','P15','P16','P17','P18', ...
    'P20','P21','P22','P23','P24','P25','P26','P27','P28','P29'}; 
subs = {'subA', 'subB'};  % subjects in each pair
runs = {'R1st','R2nd'}; % runs

%% Main
for ps = 1:p_end % pairs
    pair = pairs{ps}; % define pair ID
    for ss = 1:length(subs) % subjects
        subfolder = fullfile('..', pair, subs{ss}); % define subject folder
        load(fullfile('..', 'SPM_temp', 'Smooth_temp01.mat')) % load template matfile
        
        sum_target = []; % making empty variable for storing filenames

        for rn = 1:length(runs) % runs
            target_folder = fullfile(subfolder, runs{rn}); % define target folder
            target_listing = dir(fullfile(target_folder, 'sub*.nii'));

            % Storing filenames
            for k = 1:length(target_listing)
                target_f = [fullfile(target_folder, target_listing(k).name) ',1'];
                sum_target = [sum_target; target_f];
            end
        end % runs

        % Store variables into a structure
        matlabbatch{1, 1}.spm.spatial.smooth.data = cellstr(sum_target); % for smooth
        matlabbatch{1, 1}.spm.spatial.smooth.fwhm = [FWHMsize FWHMsize FWHMsize]; % FWHM
        matlabbatch{1, 1}.spm.spatial.smooth.prefix = ['s' num2str(FWHMsize)];

        aDate = char(datetime('now', 'Format', 'yyyyMMdd'));
        bname = fullfile('..', pair,  subs{ss}, 'batch', ['Smooth_' aDate '.mat']); % batch filename for saving
        save(bname, 'matlabbatch') % save structure
        spm_jobman('run', matlabbatch); % Run matlabbatch
        
    end % subjects
    msg = sprintf('Data of %s was processed', pair); disp(msg); % pring message on command window
end % pairs