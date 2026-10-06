%RUN_CHECKS  Reliability verification battery for the estimation code.
%
% Goes beyond the  data checks: verifies that the estimators behave
% reliably, not just that they produce good numbers on the nominal run.
%   1. Determinism (no hidden state/randomness)
%   2. NLS Gauss-Newton convergence statistics
%   3. EKF innovation consistency (NIS of accepted measurements vs dof)
%   4. All 12 configs: estimates finite and inside the tunnel ( check 4)
%   5. EKF robustness to initialization (vy0, P0)
%   6. Sensitivity to process-noise tuning (q x1/3 .. x3 around defaults)
%   7. Covariance consistency: truth inside the 95% error ellipse
%   8. Known edge case: meas_model at rho=0 (directly under the BS)

clear; clc;
here = fileparts(mfilename('fullpath')); root = fullfile(here,'..');
addpath(fullfile(root,'src'));
D = load_data();
noise = load_noise(root);
M = metrics(); tru = D.ue_pos(1:2,:);

fprintf('=== 1. Determinism ===\n');
[ea,ia] = ekf_track(D,'cv','joint',noise);
[eb,ib] = ekf_track(D,'cv','joint',noise);
assert(isequal(ea,eb) && isequal(ia.d2,ib.d2), 'EKF not deterministic');
[na] = nls_position(D,'joint',noise); [nb] = nls_position(D,'joint',noise);
assert(isequal(na,nb), 'NLS not deterministic');
fprintf('  EKF and NLS bit-identical across runs: OK\n');

fprintf('\n=== 2. NLS convergence (maxit=30, tol=1e-6) ===\n');
for mode = {'toa','aoa','joint'}
    [~,info] = nls_position(D,mode{1},noise);
    fprintf('  %-6s iters: median %2d  max %2d  steps hitting maxit: %d/%d\n', ...
        mode{1}, median(info.iters), max(info.iters), sum(info.iters>=30), D.N);
end

fprintf('\n=== 3. EKF innovation consistency (mean NIS of accepted meas) ===\n');
fprintf('  (mean d2 ~ dof => R,Q consistent; << dof => conservative filter)\n');
dofs = struct('toa',1,'aoa',2,'joint',3);
for model = {'cp','cv','ca'}
    for mode = {'toa','aoa','joint'}
        [~,I] = ekf_track(D,model{1},mode{1},noise);
        d2 = I.d2(~I.gated & isfinite(I.d2));
        fprintf('  %-3s %-6s mean NIS %6.2f (dof %d), gated %2d/%d\n', ...
            model{1}, mode{1}, mean(d2), dofs.(mode{1}), sum(I.gated), D.N);
    end
end

fprintf('\n=== 4. All 12 configs: finite + inside tunnel (Sec.7 check 4) ===\n');
R = run_estimators(D, noise, 'all');
for i = 1:numel(R)
    est = R(i).est;
    assert(all(isfinite(est(:))), 'check4:%s has NaN/Inf', R(i).label);
    assert(all(est(1,:) >= -5  & est(1,:) <= 5), 'check4:%s x out of tunnel', R(i).label);
    assert(all(est(2,:) >= -10 & est(2,:) <= 105),'check4:%s y out of tunnel', R(i).label);
end
fprintf('  all 12 estimates finite and within bounds: OK\n');

fprintf('\n=== 5. EKF (cv,joint) robustness to initialization ===\n');
for vy0 = [0 13 27 40]
    est = ekf_track(D,'cv','joint',noise,struct('vy0',vy0));
    s = M.stats(M.pos_err(est,tru));
    fprintf('  vy0 = %2d m/s : RMSE %.3f m\n', vy0, s.rmse);
end
for P0p = [1 5 20]
    est = ekf_track(D,'cv','joint',noise,struct('P0pos',P0p));
    s = M.stats(M.pos_err(est,tru));
    fprintf('  P0pos = %2d m : RMSE %.3f m\n', P0p, s.rmse);
end

fprintf('\n=== 6. Sensitivity to q around tuned defaults ===\n');
cfgs = {{'cv','joint'},{'cv','aoa'},{'ca','joint'},{'cv','toa'}};
for c = cfgs
    md = c{1}{1}; mo = c{1}{2};
    q0 = ekf_default_q(md,mo);
    line = sprintf('  %-3s %-6s (q0=%g):', md, mo, q0);
    for f = [1/3 1 3]
        est = ekf_track(D,md,mo,noise,struct('q',q0*f));
        s = M.stats(M.pos_err(est,tru));
        line = [line sprintf('  x%-4.2g->%6.3f', f, s.rmse)]; %#ok<AGROW>
    end
    fprintf('%s m\n', line);
end

fprintf('\n=== 7. Covariance consistency: truth inside the 95%% error ellipse ===\n');
fprintf('  (~95%% = well-calibrated P; far below = filter overconfident)\n');
for mode = {'toa','aoa','joint'}
    [est, I] = ekf_track(D,'cv',mode{1},noise);
    inside = false(1,D.N);
    for t = 1:D.N
        P = I.Ppos(:,:,t);
        [~,~,s95] = error_ellipse(P);
        d = est(:,t) - tru(:,t);
        inside(t) = (d' * (P \ d)) <= s95;
    end
    fprintf('  cv %-6s : truth inside 95%% ellipse at %3.0f%% of steps\n', ...
        mode{1}, 100*mean(inside));
end
fprintf(['  NOTE: toa-only is far below 95%% by design - its lateral error is a\n' ...
         '        DRIFT (a bias), and a Kalman covariance only models zero-mean\n' ...
         '        random error, so P cannot represent it. The filter is not just\n' ...
         '        imprecise there, it is inconsistent. aoa/joint are ~85-90%%.\n']);

fprintf('\n=== 8. Edge case: meas_model directly above/below BS (rho=0) ===\n');
[z,H] = meas_model([D.bs(1); D.bs(2); 1.5], D.bs, -90, -30, 'joint', true);
assert(all(isfinite(z)) && all(isfinite(H(:))), 'check7:rho0', ...
    'meas_model must stay finite at rho=0 (guarded Jacobian)');
fprintf('  z and H finite at rho=0 (guarded Jacobian): OK\n');

fprintf('\nAll reliability checks completed.\n');
