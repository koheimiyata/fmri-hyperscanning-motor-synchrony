%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   B2_HRFConvolution
%
%   2025/12/06
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Convolve behavioral regressors with the hemodynamic response
%            function (HRF) for use as SPM task regressors in G1_FLevel_GLM.
%            Processing steps:
%              1. Downsample behavioral signals from 30 Hz to 1 Hz
%              2. Convolve with SPM canonical HRF (causal convolution)
%              3. Remove common HRF bias (mean across all subjects/sessions)
%                 via regression
%              4. Trim to 4-367 s and z-score normalize
%            Regressors: theta velocity, theta acceleration, approach
%            absolute acceleration, MEA, relative phase.
%
%   Note: Causal convolution is implemented by taking the first T points
%         of conv(x, hrf) rather than using conv(..., 'same').
%         The common HRF bias (scan-onset-locked component) is estimated
%         data-driven and regressed out across all subjects and sessions.
%
%   Input : FingerTip_theta.mat (FingerTip struct),
%           MotionEnergy.mat (MEA struct)
%   Output: HRF_Convolution.mat (HRF_Convolution struct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% Initial settings
pairs = {'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', ...
         'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18', ...
         'P20', 'P21', 'P22', 'P23', 'P24', 'P25', 'P26', 'P27', 'P28', 'P29'};
subjects = {'subA', 'subB'};
sessions = {'R1st', 'R2nd'};

fps = 30;   % frame rate of behavioral data (Hz)
TR  = 1;    % fMRI TR (s)
st  = 1;    % start frame (no volumes discarded at this stage)
ed  = st + fps*369 - 1;  % end frame (369 s at 30 Hz)
T   = 369;  % time series length after downsampling to 1 Hz

%% Load data
load(fullfile('..', 'Mat_file', 'FingerTip_theta.mat'))
load(fullfile('..', 'Mat_file', 'MotionEnergy.mat'))

hrf = spm_hrf(TR);

%% Preallocate accumulators for common HRF bias estimation
sum_vel = [];
sum_acc = [];
sum_abs = [];
sum_rel = [];
sum_app = [];
sum_mea = [];

HRF_raw = struct();

%% Pass 1: HRF convolution and accumulation for common bias estimation
for ps = 1:length(pairs)
    for sn = 1:2

        % ---- 1. Individual phase (Theta) and MEA ----
        for ss = 1:length(subjects)

            % Load theta and MEA at 30 Hz
            theta_30Hz = FingerTip.(pairs{ps}).(subjects{ss}).(sessions{sn}).theta;
            mea_30Hz   = MEA.(pairs{ps}).(subjects{ss}).(sessions{sn}).RAW;

            theta_30Hz = theta_30Hz(st:ed);
            mea_30Hz   = mea_30Hz(st:ed);

            % Store theta per subject for relative phase computation
            if ss == 1
                theta_A = theta_30Hz;
            else
                theta_B = theta_30Hz;
            end

            % Smooth and differentiate theta
            uw_theta30Hz = unwrap(theta_30Hz);
            theta_smooth = smooth(uw_theta30Hz, 15);
            mea_smooth   = smooth(mea_30Hz,     15);

            pre_v_theta_30Hz = diff(theta_smooth);
            v_theta_30Hz     = [pre_v_theta_30Hz(1); pre_v_theta_30Hz];  % velocity
            av_theta_30Hz    = [0; diff(v_theta_30Hz)];                  % acceleration
            aav_theta_30Hz   = abs(av_theta_30Hz);                       % absolute acceleration

            % Store acceleration per subject for approach computation
            if ss == 1
                acc_A = av_theta_30Hz(st:ed);
            else
                acc_B = av_theta_30Hz(st:ed);
            end

            % Downsample 30 Hz -> 1 Hz
            n = floor(length(v_theta_30Hz) / fps);
            v_theta_1Hz   = zeros(n, 1);
            av_theta_1Hz  = zeros(n, 1);
            aav_theta_1Hz = zeros(n, 1);
            mea_1Hz       = zeros(n, 1);

            for i = 1:n
                idx_start = (i-1)*fps + 1;
                idx_end   = i*fps;
                v_theta_1Hz(i)   = mean(v_theta_30Hz(idx_start:idx_end))   * fps;
                av_theta_1Hz(i)  = mean(av_theta_30Hz(idx_start:idx_end))  * fps;
                aav_theta_1Hz(i) = mean(aav_theta_30Hz(idx_start:idx_end)) * fps;
                mea_1Hz(i)       = mean(mea_smooth(idx_start:idx_end))     * fps;
            end

            % Detrend
            v_theta_1Hz   = detrend(v_theta_1Hz);
            av_theta_1Hz  = detrend(av_theta_1Hz);
            aav_theta_1Hz = detrend(aav_theta_1Hz);
            mea_1Hz       = detrend(mea_1Hz);

            % Causal HRF convolution (take first T points)
            reg_conv_hrf_v   = conv(v_theta_1Hz,   hrf);
            reg_conv_hrf_a   = conv(av_theta_1Hz,  hrf);
            reg_conv_hrf_aa  = conv(aav_theta_1Hz, hrf);
            reg_conv_hrf_mea = conv(mea_1Hz,       hrf);

            % Store raw convolved signals
            HRF_raw.(pairs{ps}).(sessions{sn}).Theta_velo.(subjects{ss})    = reg_conv_hrf_v(1:T);
            HRF_raw.(pairs{ps}).(sessions{sn}).Theta_acc.(subjects{ss})     = reg_conv_hrf_a(1:T);
            HRF_raw.(pairs{ps}).(sessions{sn}).Theta_abs_acc.(subjects{ss}) = reg_conv_hrf_aa(1:T);
            HRF_raw.(pairs{ps}).(sessions{sn}).MEA.(subjects{ss})           = reg_conv_hrf_mea(1:T);

            % Accumulate for common bias estimation
            sum_vel = [sum_vel, reg_conv_hrf_v(1:T)];
            sum_acc = [sum_acc, reg_conv_hrf_a(1:T)];
            sum_abs = [sum_abs, reg_conv_hrf_aa(1:T)];
            sum_mea = [sum_mea, reg_conv_hrf_mea(1:T)];

        end % subjects

        % ---- 2. Relative phase ----
        phi_rad        = angle(exp(1i*(theta_A - theta_B)));
        ab_phi_rad     = abs(phi_rad);
        phi_rad_smooth = smooth(ab_phi_rad, 15);

        % Downsample 30 Hz -> 1 Hz
        n = floor(length(phi_rad_smooth) / fps);
        phi_rad_1Hz = zeros(n, 1);
        for i = 1:n
            idx_start = (i-1)*fps + 1;
            idx_end   = i*fps;
            phi_rad_1Hz(i) = mean(phi_rad_smooth(idx_start:idx_end)) * fps;
        end
        phi_rad_1Hz = detrend(phi_rad_1Hz);

        % Causal HRF convolution
        PR_conv_hrf = conv(phi_rad_1Hz, hrf);
        HRF_raw.(pairs{ps}).(sessions{sn}).Relative_phase = PR_conv_hrf(1:T);
        sum_rel = [sum_rel, PR_conv_hrf(1:T)];

        % ---- 3. Approach absolute acceleration ----
        % Identify moments when each subject is adjusting toward synchrony
        tolerance = deg2rad(10);

        idx_A_adjusting = ((phi_rad >  tolerance) & (acc_A < 0)) | ...
                          ((phi_rad < -tolerance) & (acc_A > 0));
        idx_B_adjusting = ((phi_rad >  tolerance) & (acc_B > 0)) | ...
                          ((phi_rad < -tolerance) & (acc_B < 0));

        adjusting_acc_A = zeros(size(acc_A));
        adjusting_acc_B = zeros(size(acc_B));
        adjusting_acc_A(idx_A_adjusting) = acc_A(idx_A_adjusting);
        adjusting_acc_B(idx_B_adjusting) = acc_B(idx_B_adjusting);

        abs_adjusting_acc_A = abs(adjusting_acc_A);
        abs_adjusting_acc_B = abs(adjusting_acc_B);

        for ss = 1:length(subjects)
            if ss == 1
                abs_adjusting_acc = abs_adjusting_acc_A;
            else
                abs_adjusting_acc = abs_adjusting_acc_B;
            end

            % Downsample 30 Hz -> 1 Hz
            n = floor(length(abs_adjusting_acc) / fps);
            abs_adjusting_acc_1Hz = zeros(n, 1);
            for i = 1:n
                idx_start = (i-1)*fps + 1;
                idx_end   = i*fps;
                abs_adjusting_acc_1Hz(i) = mean(abs_adjusting_acc(idx_start:idx_end)) * fps;
            end
            abs_adjusting_acc_1Hz = detrend(abs_adjusting_acc_1Hz);

            % Causal HRF convolution
            reg_conv_hrf_aaa = conv(abs_adjusting_acc_1Hz, hrf);
            HRF_raw.(pairs{ps}).(sessions{sn}).Theta_approach_abs_acc.(subjects{ss}) = reg_conv_hrf_aaa(1:T);
            sum_app = [sum_app, reg_conv_hrf_aaa(1:T)];
        end % subjects

    end % sessions
end % pairs

%% Compute common HRF bias (mean across all pairs/sessions/subjects)
m_vel = mean(sum_vel, 2);
m_acc = mean(sum_acc, 2);
m_abs = mean(sum_abs, 2);
m_app = mean(sum_app, 2);
m_rel = mean(sum_rel, 2);
m_mea = mean(sum_mea, 2);

% Design matrix: common HRF bias + intercept
G_vel = [m_vel, ones(T,1)];
G_acc = [m_acc, ones(T,1)];
G_abs = [m_abs, ones(T,1)];
G_app = [m_app, ones(T,1)];
G_rel = [m_rel, ones(T,1)];
G_mea = [m_mea, ones(T,1)];

%% Pass 2: Remove common HRF bias, trim to 4-367 s, and z-score normalize
for ps = 1:length(pairs)
    for sn = 1:2

        % ---- Theta regressors ----
        for ss = 1:length(subjects)
            y_vel = HRF_raw.(pairs{ps}).(sessions{sn}).Theta_velo.(subjects{ss});
            y_acc = HRF_raw.(pairs{ps}).(sessions{sn}).Theta_acc.(subjects{ss});
            y_abs = HRF_raw.(pairs{ps}).(sessions{sn}).Theta_abs_acc.(subjects{ss});

            y_vel_clean = y_vel - G_vel * (G_vel \ y_vel);
            y_acc_clean = y_acc - G_acc * (G_acc \ y_acc);
            y_abs_clean = y_abs - G_abs * (G_abs \ y_abs);

            HRF_Convolution.(pairs{ps}).(sessions{sn}).Theta_velo.(subjects{ss})    = zscore(y_vel_clean(4:367));
            HRF_Convolution.(pairs{ps}).(sessions{sn}).Theta_acc.(subjects{ss})     = zscore(y_acc_clean(4:367));
            HRF_Convolution.(pairs{ps}).(sessions{sn}).Theta_abs_acc.(subjects{ss}) = zscore(y_abs_clean(4:367));
        end

        % ---- Relative phase ----
        y_rel = HRF_raw.(pairs{ps}).(sessions{sn}).Relative_phase;
        y_rel_clean = y_rel - G_rel * (G_rel \ y_rel);
        HRF_Convolution.(pairs{ps}).(sessions{sn}).Relative_phase = zscore(y_rel_clean(4:367));

        % ---- Approach absolute acceleration ----
        for ss = 1:length(subjects)
            y_app = HRF_raw.(pairs{ps}).(sessions{sn}).Theta_approach_abs_acc.(subjects{ss});
            y_app_clean = y_app - G_app * (G_app \ y_app);
            HRF_Convolution.(pairs{ps}).(sessions{sn}).Theta_approach_abs_acc.(subjects{ss}) = zscore(y_app_clean(4:367));
        end

        % ---- MEA ----
        for ss = 1:length(subjects)
            y_mea = HRF_raw.(pairs{ps}).(sessions{sn}).MEA.(subjects{ss});
            y_mea_clean = y_mea - G_mea * (G_mea \ y_mea);
            HRF_Convolution.(pairs{ps}).(sessions{sn}).MEA.(subjects{ss}) = zscore(y_mea_clean(4:367));
        end

    end % sessions
end % pairs

save(fullfile('..', 'Mat_file', 'HRF_Convolution.mat'), 'HRF_Convolution')