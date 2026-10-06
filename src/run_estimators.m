function R = run_estimators(D, noise, which)
%RUN_ESTIMATORS  Run NLS and EKF configurations and score them vs truth.
%
%   R = run_estimators(D, noise)          runs all configs
%   R = run_estimators(D, noise, which)   'nls' | 'ekf' | 'all' (default)
%
% Returns a struct array R with one entry per configuration:
%   .kind   'NLS' | 'EKF'
%   .model  '-' (NLS) | 'cp'|'cv'|'ca'
%   .mode   'toa'|'aoa'|'joint'
%   .label  short display label, e.g. 'EKF-cv joint'
%   .est    2xN estimate [x;y]
%   .err    1xN Euclidean position error vs true
%   .rmse .median .p95 .meanerr .maxerr   scalar metrics
%   .xrmse .yrmse                          per-axis RMSE
%   .gated  number of gated steps (EKF only; 0 for NLS)
%   .info   raw estimator info struct
%
% TRUE positions are used ONLY to score (.err); the estimators see only
% measured (*_hat) values. Pure function.

    if nargin < 3 || isempty(which), which = 'all'; end
    M   = metrics();
    tru = D.ue_pos(1:2,:);
    modes  = {'toa','aoa','joint'};
    models = {'cp','cv','ca'};
    R = struct([]);

    add = @(R,e) [R, e];

    if any(strcmp(which,{'nls','all'}))
        for mo = modes
            [est, info] = nls_position(D, mo{1}, noise);
            R = add(R, score('NLS','-', mo{1}, est, info, tru, M, 0));
        end
    end
    if any(strcmp(which,{'ekf','all'}))
        for md = models
            for mo = modes
                [est, info] = ekf_track(D, md{1}, mo{1}, noise);
                R = add(R, score('EKF', md{1}, mo{1}, est, info, tru, M, sum(info.gated)));
            end
        end
    end
end

function e = score(kind, model, mode, est, info, tru, M, gated)
    % check 4 as a runtime assertion: every configuration must
    % deliver finite estimates that stay inside the tunnel.
    lbl = sprintf('%s-%s %s', kind, model, mode);
    assert(all(isfinite(est(:))), 'run_estimators:finite', ...
        '%s produced non-finite estimates', lbl);
    assert(all(est(1,:) >= -5  & est(1,:) <= 5) && ...
           all(est(2,:) >= -10 & est(2,:) <= 105), ...
        'run_estimators:bounds', '%s left the tunnel bounds', lbl);

    err = M.pos_err(est, tru);
    s   = M.stats(err);
    e.kind = kind; e.model = model; e.mode = mode;
    if strcmp(model,'-')
        e.label = sprintf('NLS %s', mode);
    else
        e.label = sprintf('EKF-%s %s', model, mode);
    end
    e.est = est; e.err = err;
    e.rmse = s.rmse; e.median = s.median; e.p95 = s.p95;
    e.meanerr = s.mae; e.maxerr = s.maxabs;
    e.xrmse = sqrt(mean((est(1,:)-tru(1,:)).^2));
    e.yrmse = sqrt(mean((est(2,:)-tru(2,:)).^2));
    e.gated = gated; e.info = info;
end
