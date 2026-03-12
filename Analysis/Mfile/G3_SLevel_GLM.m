%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   G3_SLevel
%
%   2026/3/9
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Run SPM second-level one-sample t-test for each contrast.
%            Contrast images (con_*.nii) from all subjects' first-level
%            GLM results are collected and submitted to a group-level
%            one-sample t-test in SPM.
%            Contrast names are read from the CManager template to ensure
%            consistency with first-level contrasts.
%            Special characters (>, <) in contrast names are replaced
%            with (g, l) for folder name compatibility.
%
%   Input  : con_*.nii in <pair>/<sub>/Results_GLM/,
%            CManager_temp01.mat, SLevel_temp01.mat
%   Output : SPM second-level results in Group/OneT/<conName>/,
%            batch .mat file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear %   remove all variables

%% Initial setting
p_end = 27; % the number of pairs
pairs = {'P02','P03','P04','P05','P06','P07','P08','P09', ...
    'P10','P11','P12','P13','P14','P15','P16','P17','P18', ...
    'P20','P21','P22','P23','P24','P25','P26','P27','P28','P29'}; 
subs = {'subA', 'subB'};
tempNum = 1;

%% Main
% load template matfile
load(fullfile('..', 'SPM_temp', ['CManager_temp' num2str(tempNum, '%02u') '.mat']))

conNameList = [];
cEnd = length(matlabbatch{1,1}.spm.stats.con.consess);% get condition length

for con = 1:cEnd
    conNameList = [conNameList, {matlabbatch{1,1}.spm.stats.con.consess{1,con}.tcon.name}];
end

clear matlabbatch

mkdir(fullfile('..', 'Group', 'batch'))

for c = 1:cEnd

    % contrast name
    comp = ['con_' num2str(c, '%04u')];
    conName = conNameList{c};
    for i = 1:length(conName) 
        if strcmp(conName(i), '>')
            conName(i) = 'g'; % greater than
        elseif strcmp(conName(i), '<')
            conName(i) = 'l'; % less than
        end
    end
    
    load(fullfile('..', 'SPM_temp', 'SLevel_temp01.mat')) % load template matfile
    mkdir(fullfile('..', 'Group', 'OneT', conName))
    sum_target = [];

    for ps = 1:p_end
        pair = pairs{ps};       
        for ss = 1:2
            target_folder = fullfile('..', pair, subs{ss}, 'Results_GLM');
            listing_con = dir(fullfile(target_folder, [comp '.nii']));
            for k = 1:length(listing_con)
                target_f = [fullfile(target_folder, listing_con(k).name) ',1'];
                sum_target = [sum_target; target_f];
            end
        end % ss
    end % ps

    % Store variables into a structure 
    matlabbatch{1, 1}.spm.stats.factorial_design.dir = cellstr(fullfile('..', 'Group', 'OneT', conName));
    matlabbatch{1, 1}.spm.stats.factorial_design.des.t1.scans = cellstr(sum_target); %
    matlabbatch{1, 3}.spm.stats.con.consess{1,1}.tcon.name = conName;

    % batch file name for saving
    aDate = char(datetime('now', 'Format', 'yyyyMMdd'));
    bname = fullfile('..', 'Group', 'batch', ['SLevel_' aDate '_' conName '.mat']);
    save(bname, 'matlabbatch')    % save structure
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    
end % comp