function R = check_model_reproduction(D, verbose)
%CHECK_MODEL_REPRODUCTION  Measurement-model + Jacobian validation.
%
%   R = check_model_reproduction(D)
%
% For every step, evaluate meas_model at the TRUE position with the active
% array and compare the JOINT prediction to the TRUE LOS ray (the minimum
% true range and its matching angle column). Reports residual stats per
% array. Also finite-difference checks the analytic Jacobian.
%
% This validates both the measurement model and the array-switch logic.

    if nargin < 2, verbose = true; end

    rngErr = nan(1, D.N);
    azErr  = nan(1, D.N);
    elErr  = nan(1, D.N);
    rho    = nan(1, D.N);

    for t = 1:D.N
        k   = D.ant_idx(t);
        az0 = D.ant_or(k,1);  el0 = D.ant_or(k,2);
        p   = D.ue_pos(:,t);
        rho(t) = hypot(p(1)-D.bs(1), p(2)-D.bs(2));

        % TRUE LOS ray = minimum true range; same column index in aoa
        [losR, i] = min(D.toa{t});
        losAz = D.aoa{t}(1,i);
        losEl = D.aoa{t}(2,i);

        zp = meas_model(p, D.bs, az0, el0, 'joint', false);
        rngErr(t) = zp(1) - losR;
        azErr(t)  = wrap180(zp(2) - losAz);
        elErr(t)  = wrap180(zp(3) - losEl);
    end

    % Azimuth is near-singular when the vehicle is almost under the BS
    % (rho small). Validate the model in the well-posed far field and
    % report the near-BS region separately as an observability note.
    R.rho_gate = 8;                       % m, far-field threshold
    far = rho > R.rho_gate;

    R.rng = stats(rngErr);
    R.az  = stats(azErr(far));
    R.el  = stats(elErr(far));
    R.az_all = stats(azErr);
    R.near_bs_steps = find(~far);
    R.min_rho = min(rho);
    R.az_by_array = [stats(azErr(far & D.ant_idx==1)); stats(azErr(far & D.ant_idx==2))];

    %  analytic vs finite-difference Jacobian (t=1, joint, full 3D) 
    t = 1; k = D.ant_idx(t);
    az0 = D.ant_or(k,1); el0 = D.ant_or(k,2);
    p0 = D.ue_pos(:,t);
    [~, Ha] = meas_model(p0, D.bs, az0, el0, 'joint', false);
    Hn = zeros(3,3); h = 1e-6;
    for j = 1:3
        dp = zeros(3,1); dp(j) = h;
        zp = meas_model(p0+dp, D.bs, az0, el0, 'joint', false);
        zm = meas_model(p0-dp, D.bs, az0, el0, 'joint', false);
        Hn(:,j) = (zp - zm) / (2*h);
    end
    R.jac_max_abs_err = max(abs(Ha(:) - Hn(:)));

    if verbose
        fprintf('[model] LOS reproduction (pred - true LOS), far field rho>%g m (%d/%d steps):\n', ...
            R.rho_gate, sum(far), D.N);
        fprintf('        range : bias %+.4f m   std %.4f   max|.| %.4f  (all steps)\n', R.rng.bias, R.rng.std, R.rng.maxabs);
        fprintf('        az    : bias %+.4f deg std %.4f   max|.| %.4f\n', R.az.bias,  R.az.std,  R.az.maxabs);
        fprintf('        el    : bias %+.4f deg std %.4f   max|.| %.4f\n', R.el.bias,  R.el.std,  R.el.maxabs);
        fprintf('        az by array : A1 bias %+.4f std %.4f | A2 bias %+.4f std %.4f\n', ...
            R.az_by_array(1).bias, R.az_by_array(1).std, R.az_by_array(2).bias, R.az_by_array(2).std);
        fprintf('[model] near-BS singularity: min rho=%.2f m, az ill-posed at steps %s\n', ...
            R.min_rho, mat2str(R.near_bs_steps));
        fprintf('[model] analytic vs finite-diff Jacobian max abs err = %.3e\n', R.jac_max_abs_err);
    end

    % Assertions: model must reproduce TRUE LOS in the well-posed regime.
    assert(R.rng.maxabs < 0.5,  'model:range', 'range reproduction off by %.3f m', R.rng.maxabs);
    assert(R.az.std    < 1.0,   'model:az', 'far-field azimuth std %.3f deg too large', R.az.std);
    assert(R.el.std    < 1.0,   'model:el', 'far-field elevation std %.3f deg too large', R.el.std);
    assert(R.jac_max_abs_err < 1e-4, 'model:jacobian', ...
        'analytic Jacobian disagrees with finite difference (%.3e)', R.jac_max_abs_err);
end

function s = stats(e)
    e = e(~isnan(e));
    s.bias = mean(e); s.std = std(e); s.maxabs = max(abs(e)); s.rmse = sqrt(mean(e.^2));
end
