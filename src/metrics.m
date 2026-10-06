function M = metrics()
%METRICS  Library of error-metric helpers (returns a struct of handles).
%
%   M = metrics();
%   s = M.stats(e);          % scalar error stats of vector e
%   [x,F] = M.cdf(e);        % empirical CDF of |e| (or e)
%   d = M.pos_err(est,tru);  % per-step Euclidean position error (row)
%
% s fields: n, bias, std, rmse, mae, median, p95, maxabs
%
% Pure functions; used by Task 1 (measurement residuals) and Task 4
% (positioning error).

    M.stats   = @stats;
    M.cdf     = @ecdf_abs;
    M.pos_err = @pos_err;
end

function s = stats(e)
    e = e(:); e = e(~isnan(e));
    s.n      = numel(e);
    s.bias   = mean(e);
    s.std    = std(e);
    s.rmse   = sqrt(mean(e.^2));
    s.mae    = mean(abs(e));
    s.median = median(abs(e));
    s.p95    = pctl(abs(e), 95);
    s.maxabs = max(abs(e));
end

function q = pctl(x, p)
    % linear-interpolation percentile, no toolbox dependency
    x = sort(x(:)); n = numel(x);
    if n == 0, q = NaN; return; end
    if n == 1, q = x; return; end
    r = p/100 * (n-1) + 1;          % 1-based rank
    lo = floor(r); hi = ceil(r);
    q = x(lo) + (r-lo)*(x(hi)-x(lo));
end

function [x, F] = ecdf_abs(e)
    e = abs(e(:)); e = e(~isnan(e));
    x = sort(e);
    F = (1:numel(x))' / numel(x);
end

function d = pos_err(est, tru)
    % est, tru: 2xN or 3xN; returns 1xN Euclidean error over common rows
    n = min(size(est,1), size(tru,1));
    d = vecnorm(est(1:n,:) - tru(1:n,:));
end
