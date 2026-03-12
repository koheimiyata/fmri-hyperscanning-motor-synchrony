%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   P3_SaveConfoundsAsMat
%
%   2026/02/06
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Extract confound regressors from fMRIPrep confound TSV files
%            and save them in SPM-compatible formats.
%            The first 3 timepoints are excluded to match the discarded
%            scans removed in Discard_initialFiles.
%
%           Saved outputs:
%           1. Tab-delimited .txt (367 timepoints): 24 motion parameters + global signal
%               -> used as SPM multiple regressors file for CPM
%           2. Tab-delimited .txt (364 timepoints): same as above, last 3 points removed
%               -> used as SPM multiple regressors file for GLM
%           3. .mat file (Noise struct): white matter and CSF signals (367 timepoints)
%               per pair/subject/run -> used for nuisance regression
%
%   Input  : fMRIPrep confound TSV files
%   Output : *_24motion-gs.txt and *_24motion-gs_364.txt per run, Confound_fMRIPrep.mat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%% Initial setting
p_end = 27; % the number of pairs
pairs = {'P02','P03','P04','P05','P06','P07','P08','P09', ...
    'P10','P11','P12','P13','P14','P15','P16','P17','P18', ...
    'P20','P21','P22','P23','P24','P25','P26','P27','P28','P29'}; % Pair IDs
subs = {'subA', 'subB'};  % subjects in each pair
runs = {'R1st','R2nd'}; % runs

%% Main
D1 = {'trans', 'rot'}; % motion parameters type
D2 = {'x', 'y', 'z'}; % motion parameters direction 

for ps = 1:p_end % pairs
    pair = pairs{ps}; % define pair ID
    for ss = 1:length(subs) % subjects
        subfolder = fullfile('../../', 'BIDS', 'InterSync', 'derivatives', 'fMRIPrep', ['sub-', pair(end-1:end), subs{ss}(end)]); % access to BIDS folder
        
        for rn = 1:length(runs) % runs
            % define the tsv file name
            tsvName = fullfile(subfolder, 'func', ['sub-', pair(end-1:end), subs{ss}(end), '_task-circle_run-' num2str(rn, '%02d') '_desc-confounds_timeseries.tsv']);

            % load tsv file as variable S
            T = readtable(tsvName, "FileType","text",'Delimiter', '\t');

            % storing motion parameters excepting the first three points
            reapara = zeros(367,25);
            for d1 = 1:2
                for d2 = 1:3
                    reapara(:,12*d1+4*d2-15) = T.(D1{d1} + "_" + D2{d2})(4:end);
                    reapara(:,12*d1+4*d2-14) = T.(D1{d1} + "_" + D2{d2} + "_derivative1")(4:end);
                    reapara(:,12*d1+4*d2-13) = T.(D1{d1} + "_" + D2{d2} + "_derivative1_power2")(4:end);
                    reapara(:,12*d1+4*d2-12) = T.(D1{d1} + "_" + D2{d2} + "_power2")(4:end);
                end
            end
            reapara(:,25) = T.global_signal(4:end); % global signal

            % storing other nuisance timeseries excepting the first three points
            wm = T.white_matter(4:end); % white matter
            csf = T.csf(4:end); % csf

            % save realign parameter file as txt file for CPM
            txtName4cpm = fullfile('..', pair, subs{ss}, runs{rn}, ['rp_', pair, '_', subs{ss}, '_desc-confounds_timeseries_24motion-gs.txt']);
            writematrix(reapara, txtName4cpm,'Delimiter','\t')

            % save realign parameter file as txt file for GLM
            reapara(365:367,:) = [];
            txtName4glm = strrep(txtName4cpm, '.txt', '_364.txt');
            writematrix(reapara, txtName4glm, 'Delimiter', '\t')

            % storing variables into structure
            Noise.(pair).(subs{ss}).(runs{rn}).WM = wm;
            Noise.(pair).(subs{ss}).(runs{rn}).CSF = csf;

            clear wm csf
        end % run
    end % subject   
    msg = sprintf('Data of %s was processed', pair); disp(msg);
end % pair

% Save variable as matfile
matName = fullfile('..', 'Mat_file', 'Confound_fMRIPrep');
save(matName, 'Noise')