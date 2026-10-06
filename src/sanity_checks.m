function R = sanity_checks(D, verbose)
%SANITY_CHECKS  Runtime assertions from.
%
%   R = sanity_checks(D)            runs checks, prints a summary
%   R = sanity_checks(D, false)     silent
%
% Checks (data-only here; the meas_model reproduction check lives in
% check_model_reproduction.m once meas_model.m exists):
%   (a) min(toa{t}) ~= norm(ue_pos(:,t)-bs)   (LOS = true distance)
%   (b) antenna index changes exactly once, at t=59 -> 60
%   (c) vehicle height ~ constant z (1.5 m)
%
% Returns R with the measured tolerances so callers can log them.

    if nargin < 2, verbose = true; end

    % (a) LOS = minimum TRUE range ~ true BS-UE distance 
    losErr = zeros(1, D.N);
    for t = 1:D.N
        d = norm(D.ue_pos(:,t) - D.bs);
        losErr(t) = min(D.toa{t}) - d;
    end
    R.los_max_abs = max(abs(losErr));
    R.los_bias    = mean(losErr);
    % Reflected paths are >= LOS, so min(toa) cannot be far below the true
    % distance; allow ~0.5 m (a few steps lack a clean LOS true ray).
    assert(R.los_max_abs < 0.5, 'sanity:LOS', ...
        'max |min(toa)-trueLOS| = %.3f m exceeds 0.5 m', R.los_max_abs);

    % (b) array switch occurs exactly once, at 59->60 
    sw = find(diff(D.ant_idx) ~= 0);
    R.switch_steps = sw;
    assert(isscalar(sw) && sw == 59, 'sanity:switch', ...
        'array index must change once at t=59->60; found switches at %s', mat2str(sw));
    assert(all(D.ant_idx(1:59) == 1) && all(D.ant_idx(60:end) == 2), ...
        'sanity:switchVals', 'array assignment not 1 for 1..59 / 2 for 60..113');

    % (c) constant vehicle height 
    R.z_min = min(D.ue_pos(3,:));
    R.z_max = max(D.ue_pos(3,:));
    assert(R.z_max - R.z_min < 1e-6, 'sanity:z', ...
        'vehicle height not constant: z in [%.4f, %.4f]', R.z_min, R.z_max);
    assert(abs(R.z_min - D.z_known) < 1e-6, 'sanity:zVal', ...
        'vehicle height %.4f != assumed z_known %.4f', R.z_min, D.z_known);

    if verbose
        fprintf('[sanity] LOS check : max|min(toa)-dist| = %.4f m (bias %.4f) OK\n', ...
            R.los_max_abs, R.los_bias);
        fprintf('[sanity] switch    : single switch at t=59->60 OK\n');
        fprintf('[sanity] height    : z = %.4f m constant OK\n', R.z_min);
    end
end
