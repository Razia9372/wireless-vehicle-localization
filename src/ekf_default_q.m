function q = ekf_default_q(model, mode)
%EKF_DEFAULT_Q  Tuned process-noise intensity per (motion model, meas mode).
%
% Tuned by sweep (see the q-sensitivity figure produced by
% scripts/task3_motion_models.m and the x1/3..x3 check in
% scripts/run_checks.m). Different measurement modes need
% different q because their observability differs: ToA-only is weakly
% observable at a single BS and needs a larger q to stay agile, while the
% fully-observable joint mode is happy with less. Units differ by model
% (CP: position var [m^2]; CV: accel PSD [m^2/s^3]; CA: jerk PSD [m^2/s^5]).
    % Tuned with gating on and the inflated positioning R.
    key = [lower(model) '_' lower(mode)];
    switch key
        case 'cp_toa',   q = 10;     % CP+ToA is degenerate (no motion) - reported as such
        case 'cp_aoa',   q = 10;
        case 'cp_joint', q = 100;
        case 'cv_toa',   q = 100;
        case 'cv_aoa',   q = 100;
        case 'cv_joint', q = 300;
        case 'ca_toa',   q = 100;
        case 'ca_aoa',   q = 3000;
        case 'ca_joint', q = 3000;
        otherwise,       q = 100;
    end
end
