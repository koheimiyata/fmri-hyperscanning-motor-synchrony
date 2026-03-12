%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FLevel_GLM
%
%   2026/03/07
%
%   Written by KMiyata
%   Project: InterSync
%   Purpose: Build and run SPM first-level GLM batch for GLM analysis.
%            For each subject, smoothed 3D NIfTI volumes (up to vol 367)
%            and confound regressors are assembled into an SPM fMRI model.
%            Task regressors (HRF-convolved): velocity, approach abs acc,
%            MEA (other subject), relative phase.
%            Nuisance regressors: WM, CSF (364 timepoints).
%            Motion regressors: 24 motion params + global signal (*364.txt).
%            Model settings: 128s high-pass filter, FAST autocorrelation.
%            Results saved in <pair>/<sub>/Results_GLM/.
%
%   Input  : Smoothed NIfTI volumes, *364.txt, Confound_fMRIPrep.mat,
%            HRF_Convolution.mat
%   Output : SPM.mat (first-level model), batch .mat file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial setting
runs = {'R1st','R2nd'};
p_end = 27;
subs = {'subA', 'subB'};
pairs = {'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', ...
         'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18',  ...
         'P20', 'P21', 'P22', 'P23', 'P24', 'P25', 'P26', 'P27', 'P28', 'P29'};

%% Main
direct = pwd; % get the directory info

% load template matfile
load(fullfile(fileparts(direct), 'Mat_file', 'Confound_fMRIPrep.mat'))
load(fullfile(fileparts(direct), 'Mat_file', 'HRF_Convolution.mat'))
    
for ps = 1:p_end % pairs
    pair = pairs{ps}; % define pair ID
    pairfolder = fullfile(fileparts(direct), pair); % define pair folder
    
    for ss = 1:length(subs) % subjects
        load(fullfile(fileparts(direct), 'SPM_temp', 'FLevel_temp01.mat')) 
        
        for rn = 1:length(runs) % runs
            sum_target = []; % making empty variable
            reg_name = '';

            target_folder = fullfile(pairfolder, subs{ss}, runs{rn});
            listing_s5s = dir(fullfile(target_folder, 's5s*.nii'));
            for k = 1:length(listing_s5s)
                fname = listing_s5s(k).name;
                vol_num = str2double(fname(end-6:end-4));
                if vol_num <= 367
                    target_f = [fullfile(target_folder, fname) ',1'];
                    sum_target = [sum_target; target_f];
                end
            end
            
            listing_rp = dir(fullfile(target_folder, 'rp_*364.txt'));
            if ~isempty(listing_rp)
                reg_name = fullfile(target_folder, listing_rp(1).name);
            end    
                   
            % Store variables into a structure 
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).scans = cellstr(sum_target); % 
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).multi_reg = cellstr(reg_name); %       

            % Task (Phase velocty) Regressor
            velocity = HRF_Convolution.(pairs{ps}).(runs{rn}).Theta_velo.(subs{ss});
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(1).val = velocity;
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(1).name = 'Velocity';

            acceleration = HRF_Convolution.(pairs{ps}).(runs{rn}).Theta_approach_abs_acc.(subs{ss});
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(2).val = acceleration;
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(2).name = 'Acceleration';

            mea = HRF_Convolution.(pairs{ps}).(runs{rn}).MEA.(subs{rem(ss,2)+1}); % other subject's MEA
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(3).val = mea;
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(3).name = 'MEA';

            RP = HRF_Convolution.(pairs{ps}).(runs{rn}).Relative_phase;
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(4).val = RP;
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(4).name = 'RelativePhase';

            % Nuisance regressor
            wm_reg = Noise.(pair).(subs{ss}).(runs{rn}).WM;
            csf_reg = Noise.(pair).(subs{ss}).(runs{rn}).CSF;
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(5).val = wm_reg(1:364);
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(5).name = 'WM';
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(6).val = csf_reg(1:364);
            matlabbatch{1, 1}.spm.stats.fmri_spec.sess(rn).regress(6).name = 'CSF';
        end % runs
        
        % Define the folder name
        new_cdname = fullfile(fileparts(direct), pair, subs{ss}, 'Results_GLM');
        mkdir(new_cdname)
        
        % Be careful
        deleteName = fullfile(new_cdname, 'SPM.mat');
        delete(deleteName)
        
        % Store variables into a structure 
        matlabbatch{1, 1}.spm.stats.fmri_spec.dir = cellstr(new_cdname); % for stats(#1)
        matlabbatch{1, 1}.spm.stats.fmri_spec.cvi = 'FAST';
        matlabbatch{1, 2}.spm.stats.fmri_est.write_residuals = 1;

        % batch filename for saving
        aDate = char(datetime('now', 'Format', 'yyyyMMdd'));
        bname = fullfile(pairfolder, subs{ss}, 'batch', ['Flevel_GLM_' aDate '.mat']); 
        save(bname, 'matlabbatch') % save structure
        spm_jobman('run',matlabbatch); % Run matlabbatch

    end % subjects
    msg = sprintf('Data of %s was processed', pair); disp(msg); % pring message on command window
end % pairs