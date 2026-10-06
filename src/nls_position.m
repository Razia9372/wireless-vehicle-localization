function [est, info] = nls_position(D, mode, noise, opts)
%NLS_POSITION  Snapshot Gauss-Newton positioning.
%
%   [est, info] = nls_position(D, mode, noise, opts)
%
% Per time step: select the LOS (minimum-range) measured ray, build a
% mode-honest initial guess, and refine the horizontal position (x,y) with
% a prior-regularized Gauss-Newton iteration. z is fixed at D.z_known.
%
% Inputs:
%   D     - dataset struct from load_data
%   mode  - 'toa' | 'aoa' | 'joint'
%   noise - struct from Task 1 (sigma_r, sigma_az, sigma_el, R_*)
%   opts  - optional: .sigma_prior (m, default 50), .maxit (default 30),
%                     .tol (default 1e-6)
%
% Outputs:
%   est   - 2xN estimated [x;y] per step
%   info  - struct: .iters (1xN), .selIdx (1xN selected ray), .resid (1xN)
%
% Observability note: at a single BS, ToA-only resolves only the radial
% (range) direction; the lateral d.o.f. falls back to the prior (previous
% estimate). AoA (az+el) with known z is fully observable; joint is best.
%
% Pure function (no globals, no plotting).

    if nargin < 4, opts = struct; end
    sigma_prior = getf(opts,'sigma_prior',50);
    maxit       = getf(opts,'maxit',30);
    tol         = getf(opts,'tol',1e-6);
    % tunnel bounds (physical prior knowledge; keeps estimates in the tunnel)
    bounds      = getf(opts,'bounds',[-5 5; -10 105]);

    W       = inv(meas_cov(mode, noise));
    Wp      = (1/sigma_prior^2) * eye(2);
    angRows = angle_rows(mode);

    est    = nan(2, D.N);
    info.iters  = zeros(1,D.N);
    info.selIdx = zeros(1,D.N);
    info.resid  = nan(1,D.N);

    prev = [];                                   % previous estimate [x;y]
    for t = 1:D.N
        k = D.ant_idx(t); az0 = D.ant_or(k,1); el0 = D.ant_or(k,2);
        [z, idx] = select_measurement(D.toa_hat{t}, D.aoa_hat{t}, mode, struct('method','los'));
        info.selIdx(t) = idx;

        p0 = init_guess(mode, D, t, idx, k, az0, el0, prev);
        p_prior = p0(1:2);
        p = p0;

        for it = 1:maxit
            [zp, H] = meas_model(p, D.bs, az0, el0, mode, true);
            e = z - zp; e(angRows) = wrap180(e(angRows));
            grad = H'*W*e - Wp*(p(1:2) - p_prior);
            Hess = H'*W*H + Wp;
            dp   = Hess \ grad;
            p(1:2) = p(1:2) + dp;
            p(1) = min(max(p(1), bounds(1,1)), bounds(1,2));   % clip to tunnel
            p(2) = min(max(p(2), bounds(2,1)), bounds(2,2));
            p(3)   = D.z_known;
            info.iters(t) = it;
            if norm(dp) < tol, break; end
        end

        est(:,t)     = p(1:2);
        info.resid(t)= norm(e);
        prev         = p(1:2);
    end
end

%%%%%%
function p0 = init_guess(mode, D, t, idx, k, az0, el0, prev)
    r_los  = D.toa_hat{t}(idx);
    az_h   = D.aoa_hat{t}(1,idx);
    el_h   = D.aoa_hat{t}(2,idx);
    gAz    = az_h + az0;                 % global azimuth (deg)
    gEl    = el_h + el0;                 % global elevation (deg)
    dz     = D.z_known - D.bs(3);        % known vertical offset (<0)

    switch lower(mode)
        case 'joint'
            r  = r_los;
            p0 = D.bs + r*[cosd(gEl)*cosd(gAz); cosd(gEl)*sind(gAz); sind(gEl)];
        case 'aoa'
            % range from known height and elevation (no ToA used)
            r  = dz / sind(gEl);
            if ~isfinite(r) || r <= 0, r = r_los; end
            p0 = D.bs + r*[cosd(gEl)*cosd(gAz); cosd(gEl)*sind(gAz); sind(gEl)];
        case 'toa'
            % Single-range gives |y-50| only; the sign is ambiguous. Resolve
            % it with the KNOWN active array (array 1 covers y<50, array 2
            % y>50). Lateral x (unobservable far-field) carries from history.
            x0  = 0; if ~isempty(prev), x0 = prev(1); end
            rho = sqrt(max(r_los^2 - dz^2 - x0^2, 0));   % horizontal along-y range
            sgn = (k==1)*(-1) + (k==2)*(+1);
            p0  = [x0; D.bs(2) + sgn*rho; D.z_known];
        otherwise
            error('nls_position:mode','mode must be toa|aoa|joint');
    end
    p0(3) = D.z_known;
end

function r = angle_rows(mode)
    switch lower(mode)
        case 'toa',   r = [];
        case 'aoa',   r = [1 2];
        case 'joint', r = [2 3];
    end
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
