%%%%%%%%%%%%%%%%%%%%%%
%
%   B1_RelativePhase
%
%   2025/06/19
%   Written by KMiyata
%
%   Project: InterSync
%   Purpose: Analyze relative phase angle between paired individuals
%            during the circle-drawing task. Computes mean relative phase,
%            circular SD, and Synchronization Index (vector strength) for
%            real pairs and all pseudo pairs.
%            Analysis window: frames 91 to 91+30*364-1 (364 s at 30 Hz).
%
%   Input : FingerTip_theta.mat (FingerTip struct)
%   Output: SynchronyMetrics.mat (Synchrony struct)
%
%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% Initial settings
pairs = {'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', ...
         'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18',  ...
         'P20', 'P21', 'P22', 'P23', 'P24', 'P25', 'P26', 'P27', 'P28', 'P29'};
subs = {'subA', 'subB'};
runs = {'R1st', 'R2nd'};

fps = 30; 
st = 91;
ed = st+fps*364-1; 
 
%% Main
load(fullfile('..', 'Mat_file', 'FingerTip_theta.mat'))

rpre=[];
rpost=[];
ppre=[];
ppost=[];
F = [];
D = [];

for rn = 1:2 
    for psa = 1:length(pairs)
        % load theta of sub A
        theta_A = FingerTip.(pairs{psa}).(subs{1}).(runs{rn}).theta(st:ed);
        f_mean_A = FingerTip.(pairs{psa}).(subs{1}).(runs{rn}).mean_frequency;
        
        for psb = 1:length(pairs)
            % load theta of subB
            theta_B = FingerTip.(pairs{psb}).(subs{2}).(runs{rn}).theta(st:ed);
            f_mean_B = FingerTip.(pairs{psb}).(subs{2}).(runs{rn}).mean_frequency;
        
            phi_rad = angle(exp(1i*(theta_A - theta_B)));
            
            vec_x = cos(phi_rad);
            vec_y = sin(phi_rad);
            m_vec_x = mean(vec_x);
            m_vec_y = mean(vec_y);
            m_r_deg=atan2d(m_vec_y, m_vec_x);
    
            % calculate SD
            rel_r = sqrt(m_vec_x^2+m_vec_y^2);
            s_rp = (180/pi)*sqrt(2*(1-rel_r));
    
            % store mean and SD into variables
            if psa == psb
                rM(psa,rn) = m_r_deg;
                rSD(psa,rn) = s_rp;
                rSyncIndex(psa,rn) = rel_r;
                if rn == 1
                    rpre = [rpre; phi_rad];
                else
                    rpost = [rpost; phi_rad];
                end
                rmindex(psa) = 27*psa+psb-27;
                F = [F;f_mean_A;f_mean_B];
                D = [D;abs(f_mean_A-f_mean_B)];
            else
                pM(27*psa+psb-27,rn) = m_r_deg;
                pSD(27*psa+psb-27,rn) = s_rp;
                pSyncIndex(27*psa+psb-27,rn) = rel_r;
                if rn == 1
                    ppre = [ppre; phi_rad];
                else
                    ppost = [ppost; phi_rad];
                end
            end
        end % psb
    end % psa
end % runs

% save results as mat file
Synchrony.Real.Mean = abs(rM);
Synchrony.Real.SD = rSD;
Synchrony.Real.SyncIndex = rSyncIndex;

pM(rmindex(1:end-1),:) = [];
pSD(rmindex(1:end-1),:) = [];
pSyncIndex(rmindex(1:end-1),:) = [];

Synchrony.Pseudo.Mean = abs(pM);
Synchrony.Pseudo.SD = pSD;
Synchrony.Pseudo.SyncIndex = pSyncIndex;
save(fullfile('..', 'Mat_file', 'SynchronyMetrics.mat'), 'Synchrony')

% plot histogram
figure
width = 29.7 * 2/3;   % cm
height = 7;          % cm

set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [5 5 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);

subplot(1,2,1)
h1 = histogram(abs(rad2deg(rpre)), 18, 'Normalization','probability','FaceColor', [0.97,0.46,0.43],'EdgeColor', [0.97,0.46,0.43],'EdgeAlpha',1); hold on
h2 = histogram(abs(rad2deg(rpost)), 18, 'Normalization','probability','FaceColor', [0.78,0.49,1.00],'EdgeColor', [0.78,0.49,1.00],'EdgeAlpha',1);
ylim([0 0.2])
xticks(0:30:180);
xticklabels(string(0:30:180));
set(gca, 'FontSize', 9);
ylabel('Probability', 'FontSize', 10)
xlabel('Absolute relative phase angles (deg.)', 'FontSize', 10)
legend([h1 h2], {'Real 1st', 'Real 2nd'}, 'FontSize', 10, 'Location', 'northeast');
legend boxoff;
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none', 'XColor', 'w', 'YColor', 'w', 'LineWidth', 2);
box off;

subplot(1,2,2)
h3 = histogram(abs(rad2deg(ppre)), 18, 'Normalization','probability','FaceColor', [0.00,0.75,0.77],'EdgeColor', [0.00,0.75,0.77],'EdgeAlpha',1); hold on
h4 = histogram(abs(rad2deg(ppost)), 18, 'Normalization','probability','FaceColor', [0.49,0.68,0.00],'EdgeColor', [0.49,0.68,0.00],'EdgeAlpha',1);
ylim([0 0.2])
xticks(0:30:180);
xticklabels(string(0:30:180));
set(gca, 'FontSize', 9);
ylabel('Probability', 'FontSize', 10)
xlabel('Absolute relative phase angles (deg.)', 'FontSize', 10)
legend([h3 h4], {'Pseudo 1st', 'Pseudo 2nd'}, 'FontSize', 10, 'Location', 'northeast');
legend boxoff;
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none', 'XColor', 'w', 'YColor', 'w', 'LineWidth', 1);
box off;

%exportgraphics(gcf, 'figure.png', 'BackgroundColor', 'k', 'Resolution', 300);