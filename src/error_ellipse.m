function [ex, ey, s] = error_ellipse(P, mu, conf, npts)
%ERROR_ELLIPSE  Points of a 2-D confidence ellipse from a covariance matrix.
%
%   [ex, ey, s] = error_ellipse(P, mu, conf, npts)
%
% Inputs:
%   P    - 2x2 position covariance (m^2)
%   mu   - 2x1 centre [x; y] (default [0;0])
%   conf - confidence level (default 0.95)
%   npts - points on the ellipse (default 60)
%
% Outputs:
%   ex, ey - 1xnpts ellipse coordinates (closed curve)
%   s      - the chi-square scale used, so callers can reuse it for the
%            containment test  d' * inv(P) * d <= s
%
% The ellipse is {p : (p-mu)' P^-1 (p-mu) = s}. For 2 d.o.f. the chi-square
% CDF is 1 - exp(-x/2), so the quantile has the closed form s = -2*ln(1-conf)
% (5.991 at 95%) - no Statistics Toolbox needed.
%
% Pure function.

    if nargin < 2 || isempty(mu),   mu   = [0; 0]; end
    if nargin < 3 || isempty(conf), conf = 0.95;   end
    if nargin < 4 || isempty(npts), npts = 60;     end

    s = -2 * log(1 - conf);

    [V, Dg] = eig((P + P')/2);            % symmetrise: P is a covariance
    d = max(diag(Dg), 0);                 % clip tiny negative eigenvalues

    t    = linspace(0, 2*pi, npts);
    pts  = V * diag(sqrt(s * d)) * [cos(t); sin(t)];

    ex = pts(1,:) + mu(1);
    ey = pts(2,:) + mu(2);
end
