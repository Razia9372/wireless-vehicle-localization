%TASK4_EVALUATION  Positioning-error evaluation.
%
% Runs every configuration {NLS, EKF-CP/CV/CA} x {toa, aoa, joint}, builds a
% comparison table, CDFs of position error, an RMSE bar chart, and the hero
% est-vs-true trajectory overlay for the best config. Writes the table to
% results/metrics_table.csv.
%
% Figures (results/figures/):
%   task4_rmse_bar.png     RMSE of all configs (grouped by mode)
%   task4_cdf.png          CDF of position error, key configs
%   task4_overlay_best.png best config (EKF-CV joint) vs true
%   task4_ellipses.png     EKF-CV 95% covariance ellipses per mode (observability
%                          + filter consistency; drawn from P, no truth used)

clear; clc; close all;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..');
addpath(fullfile(root,'src'));
use_light_figures();
figdir = fullfile(root,'results','figures'); if ~isfolder(figdir), mkdir(figdir); end

D = load_data();
noise = load_noise(root);
M = metrics(); tru = D.ue_pos(1:2,:);

R = run_estimators(D, noise, 'all');     % 12 configs

% comparison table (print + CSV)
fprintf('\n=== Task 4: all configurations ===\n');
fprintf('%-16s | %7s %7s %7s %7s %7s %7s %5s\n', ...
    'config','RMSE','median','p95','mean','max','y/xRMSE','gate');
fid = fopen(fullfile(root,'results','metrics_table.csv'),'w');
fprintf(fid,'config,kind,model,mode,RMSE,median,p95,mean,max,xRMSE,yRMSE,gated\n');
for i = 1:numel(R)
    r = R(i);
    fprintf('%-16s | %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f/%4.2f %5d\n', ...
        r.label, r.rmse, r.median, r.p95, r.meanerr, r.maxerr, r.yrmse, r.xrmse, r.gated);
    fprintf(fid,'%s,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n', ...
        r.label,r.kind,r.model,r.mode,r.rmse,r.median,r.p95,r.meanerr,r.maxerr,r.xrmse,r.yrmse,r.gated);
end
fclose(fid);
fprintf('Wrote results/metrics_table.csv\n');

% index helper
findcfg = @(lab) R(find(strcmp({R.label},lab),1));

% RMSE bar chart, grouped by mode
modes = {'toa','aoa','joint'};
groups = {'NLS','EKF-cp','EKF-cv','EKF-ca'};
RM = nan(numel(modes), numel(groups));
for a = 1:numel(modes)
    RM(a,1) = findcfg(sprintf('NLS %s',modes{a})).rmse;
    RM(a,2) = findcfg(sprintf('EKF-cp %s',modes{a})).rmse;
    RM(a,3) = findcfg(sprintf('EKF-cv %s',modes{a})).rmse;
    RM(a,4) = findcfg(sprintf('EKF-ca %s',modes{a})).rmse;
end
f = figure('Color','w','Position',[100 100 900 460]);
b = bar(RM); grid on; set(gca,'YScale','log','XTickLabel',upper(modes),'FontSize',11);
ylim([0.2 150]);
ylabel('position RMSE (m, log)'); xlabel('measurement mode');
title('Task 4 - position RMSE: NLS vs EKF (CP/CV/CA)');
legend(groups,'Location','northeast');
for k=1:numel(b), xtips=b(k).XEndPoints; ytips=b(k).YEndPoints;
    text(xtips,ytips,compose('%.2f',RM(:,k)),'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','FontSize',8); end
exportgraphics(f, fullfile(figdir,'task4_rmse_bar.png'),'Resolution',150);

% CDFs of position error (key configs)
key = {'NLS joint','EKF-cv toa','EKF-cv aoa','EKF-cv joint'};
f = figure('Color','w','Position',[100 100 760 480]); hold on; box on; grid on;
co = lines(numel(key));
for i=1:numel(key)
    r = findcfg(key{i}); [x,Fc] = M.cdf(r.err);
    plot(x, Fc, 'LineWidth', 1.8, 'Color', co(i,:), 'DisplayName', key{i});
end
set(gca,'XScale','log'); xlim([1e-2 1e2]); ylim([0 1]);
xlabel('position error (m, log)'); ylabel('CDF');
title('Task 4 - CDF of position error'); legend('Location','southeast'); set(gca,'FontSize',11);
exportgraphics(f, fullfile(figdir,'task4_cdf.png'),'Resolution',150);

% hero overlay: best config vs true
best = findcfg('EKF-cv joint');
plot_traj_overlay(D, {best.est}, {sprintf('EKF-cv joint (RMSE %.2f m)',best.rmse)}, ...
    'Task 4 - best estimate vs true trajectory', ...
    fullfile(figdir,'task4_overlay_best.png'));

% covariance ellipses: what the filter THINKS it knows 
% The EKF reports a position covariance P at every step. Its 95% ellipse is
% drawn straight from the estimator's own output (no truth involved), so its
% SHAPE is a direct picture of observability: a long thin ellipse means the
% filter is admitting it cannot resolve that direction. We also score how
% often the truth actually falls inside the ellipse (consistency).
modes = {'toa','aoa','joint'};
step  = 8;                                  % draw an ellipse every 8th step
co    = lines(3);

f = figure('Color','w','Position',[100 100 1150 720]);
fprintf('\n=== EKF-cv covariance ellipses (95%%) ===\n');
fprintf('%-6s | %9s %9s | %-18s | %s\n','mode','major(m)','minor(m)','major axis along','truth inside 95%% ellipse');
for a = 1:numel(modes)
    mo = modes{a};
    [est, info] = ekf_track(D, 'cv', mo, noise);

    maj = nan(1,D.N); mnr = nan(1,D.N); tilt = nan(1,D.N); inside = false(1,D.N);
    for t = 1:D.N
        P = info.Ppos(:,:,t);
        [V, Dg] = eig((P+P')/2);
        [ev, ix] = sort(max(diag(Dg),0), 'descend'); V = V(:,ix);
        [~,~,s95] = error_ellipse(P);                    % chi2 scale, 2 dof
        maj(t) = sqrt(s95*ev(1)); mnr(t) = sqrt(s95*ev(2));
        tilt(t) = atan2d(abs(V(2,1)), abs(V(1,1)));      % 0 = lateral x, 90 = along y
        d = est(:,t) - tru(:,t);
        inside(t) = (d' * (P \ d)) <= s95;
    end
    axname = 'x (lateral)'; if median(tilt) > 45, axname = 'y (along tunnel)'; end
    fprintf('%-6s | %9.2f %9.2f | %-18s | %.0f%%\n', mo, median(maj), median(mnr), ...
        axname, 100*mean(inside));

    subplot(3,1,a); hold on; box on; grid on;
    plot(D.ue_pos(2,:), D.ue_pos(1,:), 'k-', 'LineWidth', 2, 'DisplayName','true');
    for t = 1:step:D.N
        % plot in (y, x) axes to match the other trajectory figures
        [ex, ey] = error_ellipse(info.Ppos(:,:,t), est(:,t), 0.95);
        hE = plot(ey, ex, '-', 'Color', [co(a,:) 0.85], 'LineWidth', 1.1);
        if t == 1, set(hE,'DisplayName','95% covariance ellipse'); else, set(hE,'HandleVisibility','off'); end
    end
    plot(est(2,:), est(1,:), '-', 'Color', co(a,:)*0.6, 'LineWidth', 1.2, 'DisplayName','estimate');
    plot(D.bs(2), 0, 'kp', 'MarkerSize',13,'MarkerFaceColor',[.5 .5 .5],'HandleVisibility','off');
    % axis equal is REQUIRED or the ellipses are geometrically distorted; the
    % vertical range is fitted per panel so the ellipses are not clipped.
    yl = max(7, 1.25*median(maj) + 2);
    axis equal; xlim([-8 106]); ylim([-yl yl]);
    ylabel('x  (lateral, m)');
    title(sprintf('%s  -  95%% ellipse semi-axes %.1f m x %.1f m (%.0f:1) along %s;  truth inside %.0f%% of steps', ...
        upper(mo), median(maj), median(mnr), median(maj)/median(mnr), axname, 100*mean(inside)));
    if a == 1, legend('Location','northwest','Orientation','horizontal'); end
    if a == 3, xlabel('y  (along tunnel, m)'); end
    set(gca,'FontSize',10);
end
sgtitle('Task 4 - what the filter knows it knows: EKF-CV position covariance (no truth used to draw the ellipses)');
exportgraphics(f, fullfile(figdir,'task4_ellipses.png'),'Resolution',150);

fprintf('\nSaved Task-4 figures to %s\n', figdir);
