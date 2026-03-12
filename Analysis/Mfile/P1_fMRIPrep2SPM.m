%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   P1_fMRIPrep2SPM
%
%   2026/02/06
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Convert preprocessed fMRI data from fMRIPrep output format
%            to SPM-compatible format.
%            Specifically, this script:
%              1. Gunzips .nii.gz anatomical images (T1w, GM/WM/CSF probmaps)
%                 into subject-specific T1 directories.
%              2. Gunzips .nii.gz functional (BOLD) images.
%              3. Converts 4D BOLD NIfTI files into 3D volumes using SPM's
%                 4D-to-3D utility (spm.util.split), saving them per run.
%              4. Saves the SPM batch file for traceability.
%
%   Input  : fMRIPrep derivatives (MNI152NLin2009cAsym space)
%   Output : 3D NIfTI volumes and SPM batch files in subject directories
%
%   Dependency: SPM25, fMRIPrep v21+
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial setting
p_end = 27; % the number of pairs
pairs = {'P02','P03','P04','P05','P06','P07','P08','P09', ...
    'P10','P11','P12','P13','P14','P15','P16','P17','P18', ...
    'P20','P21','P22','P23','P24','P25','P26','P27','P28','P29'}; % Pair IDs
subs = {'subA', 'subB'};  % subjects in each pair
runs = {'R1st','R2nd'}; % runs

%% Main
for ps = 1:p_end % pairs
    pair = pairs{ps}; % define pair ID
    for ss = 1:length(subs) % subjects
        load(fullfile('..', 'SPM_temp', 'F4D2T3D.mat')) % load mat file
        % duplicate processes in the template
        matlabbatch(1,2) = matlabbatch(1,1);

        % convert anatomical images 
        subfolder = fullfile('../..', 'BIDS', 'InterSync', 'derivatives', 'fMRIPrep', ['sub-', pair(end-1:end), subs{ss}(end)]); % define fMRIPrep output directory
        T1file = fullfile(subfolder, 'anat', ['sub-', pair(end-1:end), subs{ss}(end), '_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii']); % T1
        GMfile = fullfile(subfolder, 'anat', ['sub-', pair(end-1:end), subs{ss}(end), '_space-MNI152NLin2009cAsym_label-GM_probseg.nii']); % GM probmap
        WMfile = fullfile(subfolder, 'anat', ['sub-', pair(end-1:end), subs{ss}(end), '_space-MNI152NLin2009cAsym_label-WM_probseg.nii']); % WM probmap
        CSFfile = fullfile(subfolder, 'anat', ['sub-', pair(end-1:end), subs{ss}(end), '_space-MNI152NLin2009cAsym_label-CSF_probseg.nii']); % CSF probmap
        % gunzip the files
        gunzip([T1file, '.gz'], fullfile('..', pair, subs{ss}, 'T1'))
        gunzip([GMfile, '.gz'], fullfile('..', pair, subs{ss}, 'T1'))
        gunzip([WMfile, '.gz'], fullfile('..', pair, subs{ss}, 'T1'))
        gunzip([CSFfile, '.gz'], fullfile('..', pair, subs{ss}, 'T1'))

        for rn = 1:length(runs) % run
            % load EPI 4D images
            target_file = fullfile(subfolder, 'func', ['sub-', pair(end-1:end), subs{ss}(end), '_task-circle_run-' num2str(rn, '%02d') '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii']);
            
            % gunzip the file
            gunzip([target_file, '.gz']) 
            
            % Store images into a structure for converting 4D into 3D
            matlabbatch{1, rn}.spm.util.split.vol = cellstr(target_file);
            matlabbatch{1, rn}.spm.util.split.outdir = cellstr(fullfile('..', pair, subs{ss}, runs{rn}));
            
            % make directory
            mkdir(fullfile('..', pair, subs{ss}, runs{rn}))
        end % run
        
        % make batch folder for SPM analysis
        mkdir(fullfile('..', pair, subs{ss}, 'batch'))
    
        % save the batch file
        aDate = char(datetime('now', 'Format', 'yyyyMMdd'));
        bname = fullfile('..', pair, subs{ss}, 'batch', ['fMRIPrep2SPM_', aDate, '.mat']); % batch filename for saving
        save(bname, 'matlabbatch') % save structure   
        spm_jobman('run', matlabbatch); % run matlabbatch

    end % subject     
    msg = sprintf('Data of %s was processed', pair); disp(msg);
end % pair