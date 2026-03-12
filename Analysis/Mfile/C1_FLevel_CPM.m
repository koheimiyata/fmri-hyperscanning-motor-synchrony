%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C1_FLevel_CPM
%
%   2026/02/06
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Build and run SPM first-level GLM batch for CPM analysis.
%            For each subject, smoothed 3D NIfTI volumes and confound
%            regressors are assembled into an SPM fMRI model specification.
%            Nuisance regressors include:
%              - 24 motion parameters + global signal (from *gs.txt file)
%              - White matter and CSF mean signals (from Confound_fMRIPrep.mat)
%            Model settings: 128s high-pass filter, no autocorrelation modeling.
%            Results are saved in <pair>/<sub>/Results_CPM/.
%
%   Input  : Smoothed NIfTI volumes, *gs.txt files, Confound_fMRIPrep.mat
%   Output : SPM.mat (first-level model), batch .mat file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial setting
p_end = 27; % the number of pairs
pairs = {'P02','P03','P04','P05','P06','P07','P08','P09', ...
    'P10','P11','P12','P13','P14','P15','P16','P17','P18', ...
    'P20','P21','P22','P23','P24','P25','P26','P27','P28','P29'}; % Pair IDs
subs = {'subA', 'subB'};  % subjects in each pair
runs = {'R1st','R2nd'}; % runs

%% Main
% load template matfile
load(fullfile('..', 'Mat_file', 'Confound_fMRIPrep.mat'))
    
for ps = 1:p_end % pairs

    pair = pairs{ps}; % define pair ID
    pairfolder = fullfile('..', pair); % define pair folder 

    for ss = 1:length(subs) % subjects

        load(fullfile('..', 'SPM_temp', 'FLevel_temp01.mat')) 
        
        for rn = 1:length(runs) % runs
            sum_target = []; % making empty variable for storing filename
            reg_name = '';
            
            % Storing filename
            target_folder = fullfile(pairfolder, subs{ss}, runs{rn}); % check the folder
            listing_sub = dir(fullfile(target_folder, 'sub*.nii'));
            for k = 1:length(listing_sub)
                target_f = [fullfile(target_folder, listing_sub(k).name) ',1'];
                sum_target = [sum_target; target_f];
            end
            
            listing_rp = dir(fullfile(target_folder, 'rp_*gs.txt'));
            if ~isempty(listing_rp)
                reg_name = fullfile(target_folder, listing_rp(1).name);
            end
            
            % Store variables into a structure 
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).scans = cellstr(sum_target); % 
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).multi_reg = cellstr(reg_name); %       

            % Regressor
            wm_reg = Noise.(pair).(subs{ss}).(runs{rn}).WM;
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(1).val = wm_reg;
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(1).name = 'WM';
            csf_reg = Noise.(pair).(subs{ss}).(runs{rn}).CSF;
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(2).val = csf_reg;
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(2).name = 'CSF';

            clear wm_reg csf_reg
        end % runs
        
        % Define the folder name
        new_cdname = fullfile('..', pair, subs{ss}, 'Results_CPM');
        mkdir(new_cdname)
                
        % Store variables into a structure 
        matlabbatch{1, 1}.spm.stats.fmri_spec.dir = cellstr(new_cdname); % for stats(#1)

        % batch filename for saving
        aDate = char(datetime('now', 'Format', 'yyyyMMdd'));
        bname = fullfile(pairfolder, subs{ss}, 'batch', ['Flevel_CPM_' aDate '.mat']); 
        save(bname, 'matlabbatch') % save structure
        spm_jobman('run', matlabbatch); % Run matlabbatch

    end % subjects
    msg = sprintf('Data of %s was processed', pair); disp(msg);
end % pairs