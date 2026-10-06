% Compares constant-position (CP), constant-velocity (CV) and
% constant-acceleration (CA) motion models in the EKF. Primary comparison
% uses the joint measurement mode; a secondary table shows the ToA-only
% case where the motion model matters most (CP cannot follow a moving
% vehicle from range alone).
%
% Figures (results/figures/):
%   task3_traj_joint.png    true vs CP/CV/CA trajectory (joint mode)
%   task3_err_time.png      position error vs time, CP/CV/CA (joint)

clear; clc; close all;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..');
addpath(fullfile(root,'src'));
use_light_figures();
figdir = fullfile(root,'results','figures'); if ~isfolder(figdir), mkdir(figdir); end

D = load_data();
noise = load_noise(root);
M = metrics(); tru = D.ue_pos(1:2,:);
models = {'cp','cv','ca'};

% run EKF for each model, joint + toa
est = struct();
for md = models
    [est.(md{1}).joint, ij] = ekf_track(D, md{1}, 'joint', noise);
    est.(md{1}).joint_g = sum(ij.gated);
    [est.(md{1}).toa, it]   = ekf_track(D, md{1}, 'toa', noise);
    est.(md{1}).toa_g = sum(it.gated);
end

% table
fprintf('\n=== Task 3: EKF motion models ===\n');
fprintf('%-4s | %-22s | %-22s\n','mdl','joint (RMSE/med/max)','toa (RMSE/med/max)');
for md = models
    sj = M.stats(M.pos_err(est.(md{1}).joint,tru));
    st = M.stats(M.pos_err(est.(md{1}).toa,tru));
    fprintf('%-4s | %6.3f %6.3f %6.3f      | %6.3f %6.3f %6.3f\n', upper(md{1}), ...
        sj.rmse,sj.median,sj.maxabs, st.rmse,st.median,st.maxabs);
end
fprintf('(process-noise q per config from ekf_default_q.m)\n');

% trajectory overlay (joint)
plot_traj_overlay(D, {est.cp.joint, est.cv.joint, est.ca.joint}, ...
    {'EKF-CP','EKF-CV','EKF-CA'}, ...
    'Task 3 - motion models (joint mode): true vs CP/CV/CA', ...
    fullfile(figdir,'task3_traj_joint.png'));

% effect of process-noise tuning (the spec asks to report this)
qs = [1 10 100 300 1000 3000 10000];
rm = nan(numel(models), numel(qs));
for i = 1:numel(models)
    for j = 1:numel(qs)
        e = ekf_track(D, models{i}, 'joint', noise, struct('q', qs(j)));
        s = M.stats(M.pos_err(e, tru));
        rm(i,j) = s.rmse;
    end
end
co = lines(3); mk = {'o-','s-','^-'};
f = figure('Color','w','Position',[100 100 760 460]); hold on; box on; grid on;
for i = 1:numel(models)
    plot(qs, rm(i,:), mk{i}, 'Color',co(i,:), 'LineWidth',1.6, ...
        'MarkerFaceColor',co(i,:), 'DisplayName', upper(models{i}));
    q0 = ekf_default_q(models{i}, 'joint');
    e  = ekf_track(D, models{i}, 'joint', noise, struct('q', q0));
    s  = M.stats(M.pos_err(e, tru));
    plot(q0, s.rmse, 'p', 'MarkerSize',15, 'Color',co(i,:), ...
        'MarkerFaceColor','w', 'LineWidth',1.4, 'HandleVisibility','off');
end
set(gca,'XScale','log','YScale','log');
xlabel('process-noise intensity q'); ylabel('position RMSE (m)');
title('Task 3 - effect of process-noise tuning (joint mode; stars = chosen q)');
legend('Location','best'); set(gca,'FontSize',11);
exportgraphics(f, fullfile(figdir,'task3_q_sweep.png'),'Resolution',150);

% error vs time (joint)
sw = 59.5; co = lines(3);
f = figure('Color','w','Position',[100 100 1000 460]); hold on; box on; grid on;
plot(1:D.N, M.pos_err(est.cp.joint,tru), '-','Color',co(1,:),'LineWidth',1.4,'DisplayName','CP');
plot(1:D.N, M.pos_err(est.cv.joint,tru), '-','Color',co(2,:),'LineWidth',1.4,'DisplayName','CV');
plot(1:D.N, M.pos_err(est.ca.joint,tru), '-','Color',co(3,:),'LineWidth',1.4,'DisplayName','CA');
xline(sw,'k--','array switch','HandleVisibility','off');
set(gca,'YScale','log'); ylim([1e-2 1e1]);
xlabel('time step t'); ylabel('position error (m)');
title('Task 3 - EKF position error vs time (joint mode)');
legend('Location','best'); set(gca,'FontSize',11);
exportgraphics(f, fullfile(figdir,'task3_err_time.png'),'Resolution',150);

fprintf('\nSaved Task-3 figures to %s\n', figdir);
