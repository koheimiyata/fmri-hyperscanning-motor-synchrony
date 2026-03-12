%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   G2_CManager
%
%   2026/03/09
%   Written by KMiyata
%   Project: InterSync
%   Purpose: Run SPM Contrast Manager for each subject's first-level GLM.
%            Loads a pre-defined contrast template (CManager_temp01.mat)
%            and applies it to each subject's SPM.mat in Results_GLM/.
%
%   Input  : SPM.mat in <pair>/<sub>/Results_GLM/,
%            CManager_temp01.mat (contrast definitions)
%   Output : Contrast images (con_*.nii, spmT_*.nii) in Results_GLM/,
%            batch .mat file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial setting
p_end = 27; % the number of pairs
pairs = {'P02','P03','P04','P05','P06','P07','P08','P09', ...
    'P10','P11','P12','P13','P14','P15','P16','P17','P18', ...
    'P20','P21','P22','P23','P24','P25','P26','P27','P28','P29'}; 
subs = {'subA', 'subB'};

%% Main
for ps = 1:p_end
    for ss = 1:2
        % load template matfile
        load(fullfile('..', 'SPM_temp', 'CManager_temp01.mat'))

        % Define the folder name
        target_file = fullfile('..', pairs{ps}, subs{ss}, 'Results_GLM', 'SPM.mat');

        % Store variables into a structure 
        matlabbatch{1, 1}.spm.stats.con.spmmat = cellstr(target_file); % for

        % batch filename for saving
        aDate = char(datetime('now', 'Format', 'yyyyMMdd'));
        bname = fullfile('..', pairs{ps}, subs{ss}, 'batch', ['CManager_' aDate '.mat']);

        % run batch
        save(bname, 'matlabbatch') % save structure
        spm_jobman('run', matlabbatch); % run matlabbatch

    end % subjects
    msg = sprintf('Data of %s was processed', pairs{ps}); disp(msg);
end % pairs