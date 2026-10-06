function [z, idx, info] = select_measurement(toa_hat, aoa_hat, mode, opts)
%SELECT_MEASUREMENT  Pick which multipath ray feeds the estimator.
%
%   [z, idx, info] = select_measurement(toa_hat, aoa_hat, mode, opts)
%
% Inputs:
%   toa_hat - 1xm measured ranges (m)
%   aoa_hat - 2xm measured angles [az;el] (deg)   (co-indexed with toa_hat)
%   mode    - 'toa' | 'aoa' | 'joint'  (which rows form z)
%   opts    - struct:
%       .method = 'los'  (default) : minimum-range measured ray (NLS)
%               = 'gate'           : Mahalanobis-nearest gated ray (EKF)
%     For 'gate':
%       .pred   - predicted measurement z_pred (mode-shaped column)
%       .Sinv   - inverse innovation covariance S^-1 (matches z size)
%       .gate   - chi-square threshold (reject if min d^2 > gate)
%
% Outputs:
%   z    - selected measurement, mode-shaped column:
%            'toa'   -> r            (1x1)
%            'aoa'   -> [az; el]     (2x1)
%            'joint' -> [r; az; el]  (3x1)
%   idx  - index of the selected ray (empty if gated out)
%   info - struct: .d2 (Mahalanobis distance of selection), .gated (logical)
%
% Pure function. Angle innovations are wrapped to [-180,180] for gating.

    if nargin < 4 || isempty(opts), opts = struct; end
    if ~isfield(opts,'method'), opts.method = 'los'; end

    m = numel(toa_hat);
    info = struct('d2', NaN, 'gated', false);

    if m == 0
        z = []; idx = []; info.gated = true; return;
    end

    switch lower(opts.method)
        case 'los'
            [~, idx] = min(toa_hat);
            z = build_z(idx);

        case 'gate'
            angRows = angle_rows(mode);
            best = inf; idx = [];
            for i = 1:m
                zc = build_z(i);
                e  = zc - opts.pred;
                e(angRows) = wrap180(e(angRows));
                d2 = e' * opts.Sinv * e;
                if d2 < best, best = d2; idx = i; end
            end
            info.d2 = best;
            if best > opts.gate
                info.gated = true; z = []; idx = [];
            else
                z = build_z(idx);
            end

        otherwise
            error('select_measurement:method', 'method must be los|gate');
    end

    %  nested: build mode-shaped measurement for ray i 
    function z = build_z(i)
        switch lower(mode)
            case 'toa',   z = toa_hat(i);
            case 'aoa',   z = aoa_hat(:,i);
            case 'joint', z = [toa_hat(i); aoa_hat(:,i)];
            otherwise, error('select_measurement:mode','mode must be toa|aoa|joint');
        end
    end
end

function r = angle_rows(mode)
    switch lower(mode)
        case 'toa',   r = [];        % range only
        case 'aoa',   r = [1 2];     % az,el
        case 'joint', r = [2 3];     % range,az,el
    end
end
