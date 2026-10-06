function [est, info] = ekf_track(D, model, mode, noise, opts)
%EKF_TRACK  Extended Kalman filter for tunnel positioning (Tasks 2 & 3).
%
%   [est, info] = ekf_track(D, model, mode, noise, opts)
%
% Inputs:
%   D     - dataset struct from load_data
%   model - motion model 'cp' | 'cv' | 'ca'
%   mode  - measurement mode 'toa' | 'aoa' | 'joint'
%   noise - noise struct from Task 1 (R_toa/R_aoa/R_joint)
%   opts  - optional:
%       .q      process-noise intensity (default per model)
%       .vy0    initial along-tunnel speed (m/s, default 27)
%       .P0pos,.P0vel,.P0acc  initial std dev of the position/velocity/
%               acceleration states (defaults 5 m, 15 m/s, 10 m/s^2)
%       .gate   chi-square gate (default 0.99 quantile for the mode dof)
%       .bounds tunnel bounds [xmin xmax; ymin ymax]
%
% Recursion: predict -> gate/select measurement -> update,
% with angle innovations wrapped to [-180,180]. Position is clipped to the
% tunnel. z is fixed at D.z_known.
%
% Outputs:
%   est  - 2xN estimated [x;y]
%   info - struct: .gated (1xN logical, true if no meas used),
%                  .d2 (1xN Mahalanobis dist of used meas),
%                  .innov (cell), .selIdx (1xN), .P (final cov), .X (full state),
%                  .Ppos (2x2xN posterior position covariance, for error_ellipse)
%
% Pure function (no globals, no plotting).

    if nargin < 5, opts = struct; end
    q      = getf(opts,'q', ekf_default_q(model, mode));
    vy0    = getf(opts,'vy0', 27);
    bounds = getf(opts,'bounds',[-5 5; -10 105]);
    gate   = getf(opts,'gate', chi2q(meas_dim(mode)));

    mm = motion_models(model, D.Ts, q);
    R  = meas_cov(mode, noise);
    angRows = angle_rows(mode);

    % initial state 
    p0 = init_position(D, 1);                  % [x;y] from first LOS ray
    x  = zeros(mm.dim,1);
    x(mm.posIdx) = p0;
    if ~isempty(mm.velIdx), x(mm.velIdx) = [0; vy0]; end   % vx0=0, vy0 along tunnel

    P0p = getf(opts,'P0pos',5)^2;
    P0v = getf(opts,'P0vel',15)^2;
    P0a = getf(opts,'P0acc',10)^2;
    Pd  = P0a*ones(mm.dim,1);
    Pd(mm.posIdx) = P0p;
    if ~isempty(mm.velIdx), Pd(mm.velIdx) = P0v; end
    P = diag(Pd);

    est = nan(2, D.N);
    info.gated  = false(1,D.N);
    info.d2     = nan(1,D.N);
    info.selIdx = nan(1,D.N);
    info.innov  = cell(1,D.N);
    info.X      = nan(mm.dim, D.N);
    info.Ppos   = nan(2, 2, D.N);

    for t = 1:D.N
        % predict 
        x = mm.F * x;
        P = mm.F * P * mm.F' + mm.Q;

        % measurement prediction at predicted position 
        k = D.ant_idx(t); az0 = D.ant_or(k,1); el0 = D.ant_or(k,2);
        p = [x(mm.posIdx); D.z_known];
        [zpred, Hp] = meas_model(p, D.bs, az0, el0, mode, true);
        H = mm.lift(Hp);
        S = H * P * H' + R;

        % gate + select 
        sel = struct('method','gate','pred',zpred,'Sinv',inv(S),'gate',gate);
        [z, idx, gi] = select_measurement(D.toa_hat{t}, D.aoa_hat{t}, mode, sel);
        info.d2(t) = gi.d2; info.selIdx(t) = iff(isempty(idx),NaN,idx);

        if isempty(z)
            info.gated(t) = true;          % no admissible measurement: predict only
        else
            % update 
            e = z - zpred; e(angRows) = wrap180(e(angRows));
            K = (P * H') / S;
            x = x + K * e;
            P = (eye(mm.dim) - K * H) * P;
            P = (P + P')/2;                % keep symmetric
            info.innov{t} = e;
        end

        % clip position to tunnel 
        x(mm.posIdx(1)) = min(max(x(mm.posIdx(1)), bounds(1,1)), bounds(1,2));
        x(mm.posIdx(2)) = min(max(x(mm.posIdx(2)), bounds(2,1)), bounds(2,2));

        est(:,t) = x(mm.posIdx);
        info.X(:,t) = x;
        info.Ppos(:,:,t) = P(mm.posIdx, mm.posIdx);   % posterior position cov
    end
    info.P = P; info.q = q; info.gate = gate; info.model = model; info.mode = mode;
end

function p0 = init_position(D, t)
    [~, i] = min(D.toa_hat{t});
    k = D.ant_idx(t); az0 = D.ant_or(k,1); el0 = D.ant_or(k,2);
    r  = D.toa_hat{t}(i);
    gAz = D.aoa_hat{t}(1,i) + az0;
    gEl = D.aoa_hat{t}(2,i) + el0;
    p = D.bs + r*[cosd(gEl)*cosd(gAz); cosd(gEl)*sind(gAz); sind(gEl)];
    p0 = p(1:2);
end

function d = meas_dim(mode)
    switch lower(mode), case 'toa', d=1; case 'aoa', d=2; case 'joint', d=3; end
end

function r = angle_rows(mode)
    switch lower(mode)
        case 'toa',   r = [];
        case 'aoa',   r = [1 2];
        case 'joint', r = [2 3];
    end
end

function g = chi2q(dof)
    % 0.99 quantile of chi-square for dof = 1..3
    tbl = [6.635, 9.210, 11.345];
    g = tbl(dof);
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function y = iff(c,a,b), if c, y=a; else, y=b; end, end
