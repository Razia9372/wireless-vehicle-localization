%TASK1_DATA_ANALYSIS  Measurement data analysis
%
% Loads data, enforces sanity + model checks, draws the scenario, then
% associates MEASURED rays to TRUE rays per step and characterises the
% measurement accuracy (range / azimuth / elevation). Derives the noise
% model R used by the NLS/EKF estimators and saves it to results/.
%
% Figures (results/figures/):
%   task1_scenario.png      tunnel + true trajectory
%   task1_raycounts.png     #rays per step, true vs measured, switch marked
%   task1_resid_hist.png    residual histograms (range/az/el)
%   task1_resid_cdf.png     empirical CDFs of |residual|
%   task1_resid_time.png    residual vs time, array switch marked

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..');
addpath(fullfile(root, 'src'));
use_light_figures();
figdir = fullfile(root, 'results', 'figures');
if ~isfolder(figdir), mkdir(figdir); end

% load + checks
D = load_data();
fprintf('Loaded %d steps, Ts=%.3f s, BS=[%g %g %g]\n', D.N, D.Ts, D.bs);
sanity_checks(D);
check_model_reproduction(D);

dist = vecnorm(diff(D.ue_pos,1,2)); spd = dist / D.Ts;
fprintf('Path %.1f m, mean speed %.1f m/s (%.0f km/h)\n', sum(dist), mean(spd), mean(spd)*3.6);

% scenario figure 
plot_scenario(D, fullfile(figdir, 'task1_scenario.png'));

% associate measured rays to true rays, collect residuals
M = metrics();
scales = [0.5; 3; 3];                 % normalising scales [m, deg, deg]
res = struct('t',[],'arr',[],'rng',[],'az',[],'el',[],'isLOS',[],'rho',[]);
ctt = zeros(1,D.N); ctm = zeros(1,D.N);

for t = 1:D.N
    k = D.ant_idx(t);
    Zt = [D.toa{t};     D.aoa{t}];        % 3xM true  [r;az;el]
    Zm = [D.toa_hat{t}; D.aoa_hat{t}];    % 3xm measured
    ctt(t) = size(Zt,2); ctm(t) = size(Zm,2);

    [pairs, ~] = associate_rays(Zt, Zm, scales);
    [~, iLOStrue] = min(D.toa{t});        % true LOS = min true range
    rho = hypot(D.ue_pos(1,t)-D.bs(1), D.ue_pos(2,t)-D.bs(2));

    for q = 1:size(pairs,1)
        im = pairs(q,1); jt = pairs(q,2);
        res.t(end+1,1)   = t;
        res.arr(end+1,1) = k;
        res.rng(end+1,1) = Zm(1,im) - Zt(1,jt);
        res.az(end+1,1)  = wrap180(Zm(2,im) - Zt(2,jt));
        res.el(end+1,1)  = wrap180(Zm(3,im) - Zt(3,jt));
        res.isLOS(end+1,1) = (jt == iLOStrue);
        res.rho(end+1,1)   = rho;
    end
end

% residual statistics 
isLOS = res.isLOS == 1;
farAz = res.rho > 8;                 % azimuth meaningful away from BS
fprintf('\n=== Measurement residuals (measured - true) ===\n');
print_stat('range  ALL', M.stats(res.rng), 'm');
print_stat('range  LOS', M.stats(res.rng(isLOS)), 'm');
print_stat('az     ALL', M.stats(res.az), 'deg');
print_stat('az LOS far', M.stats(res.az(isLOS & farAz)), 'deg');
print_stat('el     ALL', M.stats(res.el), 'deg');
print_stat('el     LOS', M.stats(res.el(isLOS)), 'deg');

% noise model R (from LOS residuals the estimator will use) 
sR   = M.stats(res.rng(isLOS));
sAz  = M.stats(res.az(isLOS & farAz));
sEl  = M.stats(res.el(isLOS));
noise.sigma_r  = sR.std;             % m
noise.sigma_az = sAz.std;            % deg
noise.sigma_el = sEl.std;            % deg
noise.R_toa    = noise.sigma_r^2;
noise.R_aoa    = diag([noise.sigma_az^2, noise.sigma_el^2]);
noise.R_joint  = diag([noise.sigma_r^2, noise.sigma_az^2, noise.sigma_el^2]);
noise.scales   = scales;
save(fullfile(root,'results','noise_model.mat'), 'noise');
fprintf('\nNoise model: sigma_r=%.3f m, sigma_az=%.3f deg, sigma_el=%.3f deg  -> saved\n', ...
    noise.sigma_r, noise.sigma_az, noise.sigma_el);

% FIGURES 
sw = 59.5;   % array switch between t=59 and 60

% (1) ray counts per step
f = figure('Color','w','Position',[100 100 900 360]); hold on; box on; grid on;
stairs(1:D.N, ctt, 'LineWidth',1.6, 'DisplayName','true rays');
stairs(1:D.N, ctm, 'LineWidth',1.6, 'DisplayName','measured rays');
xline(sw,'k--','array switch','LabelOrientation','horizontal','HandleVisibility','off');
xlabel('time step t'); ylabel('# multipath rays'); ylim([0 6]);
title('Ray count per step: true vs measured'); legend('Location','best'); set(gca,'FontSize',11);
exportgraphics(f, fullfile(figdir,'task1_raycounts.png'),'Resolution',150);

% (2) residual histograms LOS (clean) vs multipath (spread)
mp = ~isLOS;                          % non-LOS = reflected multipath rays
f = figure('Color','w','Position',[100 100 1100 340]);
hsub(1) = subplot(1,3,1);
histogram(res.rng(mp),'BinWidth',0.5,'FaceColor',[0.7 0.7 0.7],'DisplayName','multipath'); hold on;
histogram(res.rng(isLOS),'BinWidth',0.5,'FaceColor',[0 0.45 0.74],'DisplayName','LOS'); grid on;
xlim([-6 6]); xlabel('range residual (m)'); ylabel('count'); title('Range'); legend('Location','northeast');
hsub(2) = subplot(1,3,2);
histogram(res.az(mp&farAz),'BinWidth',1,'FaceColor',[0.7 0.7 0.7],'DisplayName','multipath'); hold on;
histogram(res.az(isLOS&farAz),'BinWidth',1,'FaceColor',[0.85 0.33 0.10],'DisplayName','LOS'); grid on;
xlim([-15 15]); xlabel('azimuth residual (deg)'); title('Azimuth (rho>8 m)'); legend('Location','northeast');
hsub(3) = subplot(1,3,3);
histogram(res.el(mp),'BinWidth',1,'FaceColor',[0.7 0.7 0.7],'DisplayName','multipath'); hold on;
histogram(res.el(isLOS),'BinWidth',1,'FaceColor',[0.47 0.67 0.19],'DisplayName','LOS'); grid on;
xlim([-15 15]); xlabel('elevation residual (deg)'); title('Elevation'); legend('Location','northeast');
sgtitle('Measurement residuals: LOS (min-range ray) vs multipath');
exportgraphics(f, fullfile(figdir,'task1_resid_hist.png'),'Resolution',150);

% (3) empirical CDFs of |residual| (LOS rays)
f = figure('Color','w','Position',[100 100 1100 320]);
subplot(1,3,1); [x,Fr]=M.cdf(res.rng(isLOS)); plot(x,Fr,'LineWidth',1.8); grid on;
xlabel('|range error| (m)'); ylabel('CDF'); title('Range (LOS)'); ylim([0 1]);
subplot(1,3,2); [x,Fa]=M.cdf(res.az(isLOS&farAz)); plot(x,Fa,'LineWidth',1.8); grid on;
xlabel('|az error| (deg)'); title('Azimuth (LOS, rho>8 m)'); ylim([0 1]);
subplot(1,3,3); [x,Fe]=M.cdf(res.el(isLOS)); plot(x,Fe,'LineWidth',1.8); grid on;
xlabel('|el error| (deg)'); title('Elevation (LOS)'); ylim([0 1]);
sgtitle('Empirical CDFs of LOS measurement error');
exportgraphics(f, fullfile(figdir,'task1_resid_cdf.png'),'Resolution',150);

% (4) residual vs time (LOS rays)
f = figure('Color','w','Position',[100 100 950 640]);
subplot(3,1,1); plot(res.t(isLOS),res.rng(isLOS),'.','MarkerSize',9); grid on;
xline(sw,'k--','HandleVisibility','off'); ylabel('range err (m)'); title('LOS measurement residuals vs time');
subplot(3,1,2); plot(res.t(isLOS),res.az(isLOS),'.','MarkerSize',9,'Color',[0.85 0.33 0.10]); grid on;
xline(sw,'k--','HandleVisibility','off'); ylabel('az err (deg)');
subplot(3,1,3); plot(res.t(isLOS),res.el(isLOS),'.','MarkerSize',9,'Color',[0.47 0.67 0.19]); grid on;
xline(sw,'k--','HandleVisibility','off'); ylabel('el err (deg)'); xlabel('time step t');
exportgraphics(f, fullfile(figdir,'task1_resid_time.png'),'Resolution',150);

fprintf('\nSaved Task-1 figures to %s\n', figdir);

% local helper 
function print_stat(name, s, unit)
    fprintf('  %-11s n=%3d  bias %+8.4f  std %7.4f  rmse %7.4f  p95 %7.4f  %s\n', ...
        name, s.n, s.bias, s.std, s.rmse, s.p95, unit);
end
