function [z_pred, H] = meas_model(p, bs, az0, el0, mode, fixZ)
%MEAS_MODEL  Predicted ToA/AoA measurement and Jacobian for one ray.
%
%   [z_pred, H] = meas_model(p, bs, az0, el0, mode, fixZ)
%
% Inputs:
%   p     - 3x1 vehicle position [x;y;z] (m). When fixZ is true, p(3) must
%           hold the known height (1.5 m); it is used but not estimated.
%   bs    - 3x1 base-station position (m)
%   az0   - active array boresight azimuth (deg)
%   el0   - active array boresight elevation (deg)
%   mode  - 'toa' | 'aoa' | 'joint'
%   fixZ  - logical (default true): if true, H has 2 columns (d/dx, d/dy);
%           if false, 3 columns (d/dx, d/dy, d/dz).
%
% Outputs ( VERIFIED model - "global angle minus boresight"):
%   z_pred - measurement vector:
%              'toa'   -> r
%              'aoa'   -> [az; el]      (deg, array-local frame)
%              'joint' -> [r; az; el]
%   H      - Jacobian d z_pred / d state, rows match z_pred, columns match
%            the estimated state (2 if fixZ else 3).
%
% Pure function. Angles are NOT wrapped here; wrap innovations at the caller.

    if nargin < 6 || isempty(fixZ), fixZ = true; end

    u  = p(:) - bs(:);
    dx = u(1); dy = u(2); dz = u(3);
    r    = norm(u);
    rho2 = dx^2 + dy^2;
    rho  = sqrt(rho2);
    d2r  = 180/pi;

    %  predicted measurements 
    range_pred = r;
    az_pred    = atan2d(dy, dx) - az0;
    el_pred    = asind(dz / r)  - el0;

    % full Jacobian rows w.r.t [x;y;z] 
    % Guarded denominators: az/el rows are singular at rho -> 0 (vehicle
    % exactly above/below the BS). Never reached in this dataset (min true
    % rho = 1.65 m) but keeps H finite for any query point.
    rj2 = max(rho2, 1e-9);  rj = sqrt(rj2);  rr = max(r, 1e-6);
    Hr  = [dx, dy, dz] / rr;
    Haz = d2r * [-dy/rj2, dx/rj2, 0];
    Hel = d2r * [-dx*dz/(rr^2*rj), -dy*dz/(rr^2*rj), rj/rr^2];

    switch lower(mode)
        case 'toa'
            z_pred = range_pred;          H = Hr;
        case 'aoa'
            z_pred = [az_pred; el_pred];  H = [Haz; Hel];
        case 'joint'
            z_pred = [range_pred; az_pred; el_pred];  H = [Hr; Haz; Hel];
        otherwise
            error('meas_model:mode', 'mode must be toa|aoa|joint, got "%s"', mode);
    end

    if fixZ
        H = H(:, 1:2);   % drop d/dz column (z known)
    end
end
