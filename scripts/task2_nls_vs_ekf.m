%TASK2_NLS_VS_EKF  Positioning: NLS vs EKF in toa/aoa/joint
%
% Runs NLS (snapshot Gauss-Newton) and EKF (constant-velocity motion model)
% in each of the three measurement modes, compares them quantitatively and
% qualitatively. EKF uses CV here; the motion-model comparison is Task 3.
%
% Figures (results/figures/):
%   task2_traj_<mode>.png   true vs NLS vs EKF trajectory, per mode
%   task2_err_time.png      position error vs time, all 6 configs

clear; clc; close all;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..');
addpath(fullfile(root,'src'));
use_light_figures();
figdir = fullfile(root,'results','figures'); if ~isfolder(figdir), mkdir(figdir); end

D = load_data();
noise = load_noise(root);
M = metrics(); tru = D.ue_pos(1:2,:);
modes = {'toa','aoa','joint'};

% run NLS and EKF(cv) in each mode
res = struct([]);
for i = 1:numel(modes)
    mo = modes{i};
    [eN,iN] = nls_position(D, mo, noise);
    [eE,iE] = ekf_track(D, 'cv', mo, noise);
    res(i).mode = mo; res(i).nls = eN; res(i).ekf = eE; res(i).gated = sum(iE.gated);
end

% table
fprintf('\n=== Task 2: NLS vs EKF (CV) ===\n');
fprintf('%-6s %-5s | %7s %7s %7s %7s\n','mode','est','RMSE','median','p95','max');
for i = 1:numel(modes)
    sN = M.stats(M.pos_err(res(i).nls,tru));
    sE = M.stats(M.pos_err(res(i).ekf,tru));
    fprintf('%-6s %-5s | %7.3f %7.3f %7.3f %7.3f\n', res(i).mode,'NLS',sN.rmse,sN.median,sN.p95,sN.maxabs);
    fprintf('%-6s %-5s | %7.3f %7.3f %7.3f %7.3f  (gated %d)\n', res(i).mode,'EKF',sE.rmse,sE.median,sE.p95,sE.maxabs,res(i).gated);
end

% trajectory overlays per mode
for i = 1:numel(modes)
    plot_traj_overlay(D, {res(i).nls, res(i).ekf}, ...
        {sprintf('NLS %s',res(i).mode), sprintf('EKF-cv %s',res(i).mode)}, ...
        sprintf('Task 2 - %s mode: true vs NLS vs EKF', upper(res(i).mode)), ...
        fullfile(figdir, sprintf('task2_traj_%s.png', res(i).mode)));
end

% position error vs time (all 6)
sw = 59.5;
f = figure('Color','w','Position',[100 100 1000 560]);
co = lines(3);
subplot(2,1,1); hold on; box on; grid on;
for i=1:numel(modes)
    plot(1:D.N, M.pos_err(res(i).nls,tru), '-', 'Color',co(i,:), 'LineWidth',1.3, ...
        'DisplayName',sprintf('NLS %s',res(i).mode));
end
xline(sw,'k--','HandleVisibility','off');
ylabel('error (m)'); title('NLS position error vs time'); legend('Location','best');
set(gca,'YScale','log'); ylim([1e-2 1e2]);
subplot(2,1,2); hold on; box on; grid on;
for i=1:numel(modes)
    plot(1:D.N, M.pos_err(res(i).ekf,tru), '-', 'Color',co(i,:), 'LineWidth',1.3, ...
        'DisplayName',sprintf('EKF %s',res(i).mode));
end
xline(sw,'k--','HandleVisibility','off');
xlabel('time step t'); ylabel('error (m)'); title('EKF (CV) position error vs time'); legend('Location','best');
set(gca,'YScale','log'); ylim([1e-2 1e2]);
exportgraphics(f, fullfile(figdir,'task2_err_time.png'),'Resolution',150);

fprintf('\nSaved Task-2 figures to %s\n', figdir);
