%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   P2_Discard_initialFiles
%
%   2026/02/06
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Move the first 3 volumes of each run to a 'discarded' subfolder.
%            In fMRI preprocessing, the initial scans are often excluded
%            because the MRI signal has not yet reached steady state
%            (T1 equilibration). This script isolates those volumes
%            (*_00001.nii, *_00002.nii, *_00003.nii) by moving them out of
%            the main run directory rather than deleting them, preserving
%            traceability.
%
%   Input  : 3D NIfTI volumes in each run directory (output of P1_fMRIPrep2SPM)
%   Output : First 3 volumes moved to <run>/discarded/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial settings
p_end = 27;
pairs = {'P02','P03','P04','P05','P06','P07','P08','P09', ...
    'P10','P11','P12','P13','P14','P15','P16','P17','P18', ...
    'P20','P21','P22','P23','P24','P25','P26','P27','P28','P29'};
subs = {'subA', 'subB'};
runs = {'R1st','R2nd'};

%% Main
direct = pwd; % get current directory
for ps = 1:p_end % pair
    pair = pairs{ps};
    for ss = 1:length(subs)
        pairfolder = fullfile(fileparts(direct), pair, subs{ss});
        for rn = 1:length(runs)
            direct_change = fullfile(pairfolder, runs{rn});
            cd(direct_change)

            % destination foldername
            new_foldername = fullfile(pairfolder, runs{rn}, 'discarded');
            mkdir(new_foldername)

            % find the first three scans
            target1 = dir('*_00001.nii');
            target2 = dir('*_00002.nii');
            target3 = dir('*_00003.nii');

            % check 
            if isempty(target1) || isempty(target2) || isempty(target3)
                warning('Missing files in %s', direct_change);
                continue
            end

            % file names to move
            f_target1 = fullfile(pairfolder, runs{rn}, target1.name);
            f_target2 = fullfile(pairfolder, runs{rn}, target2.name);
            f_target3 = fullfile(pairfolder, runs{rn}, target3.name);

            % move files
            movefile(f_target1, new_foldername)
            movefile(f_target2, new_foldername)
            movefile(f_target3, new_foldername)

        end % run
    end % subject
    msg = sprintf('Data of %s was processed', pair); disp(msg);  
end % pair

cd(direct) % change current directory