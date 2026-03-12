%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   C2_ResidT2F
%
%   2026/02/06
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Concatenate SPM residual images (Res_*.nii) from the
%            first-level GLM (Results_CPM) into 4D NIfTI files,
%            one per run per subject, for subsequent inter-brain
%            synchrony (CPM) analysis.
%            Residuals represent BOLD signal after nuisance regression.
%            Run 1: Res_0001–Res_0367, Run 2: Res_0368–Res_0734.
%            Output 4D files are stored in Group/InterBrain/CPM_<date>/.
%
%   Input  : Res_*.nii from <pair>/<sub>/Results_CPM/
%   Output : <pair>_<sub>_<run>_R4D.nii in Group/InterBrain/CPM_<date>/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial setting
p_end = 27; % the number of pairs
pairs = {'P02','P03','P04','P05','P06','P07','P08','P09', ...
    'P10','P11','P12','P13','P14','P15','P16','P17','P18', ...
    'P20','P21','P22','P23','P24','P25','P26','P27','P28','P29'}; % Pair IDs
subs = {'subA', 'subB'};  % subjects in each pair
runs = {'R1st','R2nd'}; % runs

%% Main
% make destination folder
aDate = char(datetime('now', 'Format', 'yyyyMMdd'));
mkdir(fullfile('..', 'Group', 'InterBrain', ['CPM_' aDate]), 'R1st')
mkdir(fullfile('..', 'Group', 'InterBrain', ['CPM_' aDate]), 'R2nd')

for ps = 1:p_end % pairs

    % load template matfile
    load(fullfile('..', 'SPM_temp', 'F3D2T4D.mat'))
    matlabbatch{1, 2} = matlabbatch{1, 1}; % duplitcate the SPM process for runs

    pair = pairs{ps}; % define pair ID
    pairfolder = fullfile('..', pair); % define pair folder   

    for ss = 1:length(subs) % subjects    

        target_folder = fullfile(pairfolder, subs{ss}, 'Results_CPM');
        
        for rn = 1:length(runs)
            sum_target = []; % making empty variable for storing filename
            for tr = 367*rn-366:367*rn
                target_f = [target_folder filesep 'Res_' num2str(tr, '%04u') '.nii,1'];
                sum_target = [sum_target; target_f];
            end % trial

            % Store variables into a structure 
            matlabbatch{1, rn}.spm.util.cat.vols = cellstr(sum_target); %             
            matlabbatch{1, rn}.spm.util.cat.name = ...
                fullfile('..', 'Group', 'InterBrain', ['CPM_' aDate], runs{rn}, [pair '_' subs{ss} '_' runs{rn} '_R4D.nii']); 
       
        end % runs        
        
        % save and run matlabbatch
        bname = fullfile(pairfolder, subs{ss}, 'batch', ['ResidT2F_CPM_' aDate '.mat']); % batch filename 
        save(bname, 'matlabbatch') % save structure
        spm_jobman('run', matlabbatch); % Run matlabbatch

    end % subjects
end % pairs