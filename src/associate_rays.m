function [pairs, info] = associate_rays(Zt, Zm, scales, gate)
%ASSOCIATE_RAYS  Greedy nearest-neighbour association of two ray sets.
%
%   [pairs, info] = associate_rays(Zt, Zm, scales, gate)
%
% Inputs:
%   Zt     - 3xM reference rays, rows [range(m); az(deg); el(deg)]
%   Zm     - 3xm query rays, same row layout
%   scales - 3x1 normalising scales [sr; saz; sel] (default [0.5;3;3])
%   gate   - max normalised squared distance to accept a pair (default Inf)
%
% Cost between ray i (Zm) and j (Zt):
%   c = (dr/sr)^2 + (wrap180(daz)/saz)^2 + (wrap180(del)/sel)^2
%
% Greedy: repeatedly take the globally smallest cost, assign that
% (query,reference) pair, remove both, until one set is exhausted or every
% remaining cost exceeds the gate. Handles M ~= m (mismatched counts).
%
% Outputs:
%   pairs - K x 2, each row [query_idx(Zm), ref_idx(Zt)]
%   info  - struct: .cost (K x1 accepted costs),
%                   .unmatched_m (query idx with no match),
%                   .unmatched_t (reference idx with no match)
%
% Pure function. No toolbox dependency.

    if nargin < 3 || isempty(scales), scales = [0.5; 3; 3]; end
    if nargin < 4 || isempty(gate),   gate   = Inf;          end

    M = size(Zt, 2);  m = size(Zm, 2);
    pairs = zeros(0, 2);  acc = [];

    % full cost matrix (m x M)
    C = inf(m, M);
    for i = 1:m
        for j = 1:M
            dr  = (Zm(1,i) - Zt(1,j)) / scales(1);
            daz = wrap180(Zm(2,i) - Zt(2,j)) / scales(2);
            del = wrap180(Zm(3,i) - Zt(3,j)) / scales(3);
            C(i,j) = dr^2 + daz^2 + del^2;
        end
    end

    rowsOpen = true(m,1);  colsOpen = true(M,1);
    while any(rowsOpen) && any(colsOpen)
        sub = C;
        sub(~rowsOpen, :) = inf;  sub(:, ~colsOpen) = inf;
        [cmin, idx] = min(sub(:));
        if ~isfinite(cmin) || cmin > gate, break; end
        [i, j] = ind2sub([m M], idx);
        pairs(end+1,:) = [i j];   %#ok<AGROW>
        acc(end+1,1)   = cmin;    %#ok<AGROW>
        rowsOpen(i) = false;  colsOpen(j) = false;
    end

    info.cost        = acc;
    info.unmatched_m = find(rowsOpen);
    info.unmatched_t = find(colsOpen);
end
