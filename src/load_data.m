function D = load_data(matfile)
%LOAD_DATA  Load the LNSM 2026 dataset into a clean struct.
%
%   D = LOAD_DATA()           loads data/LNSM_project_data_2026.mat
%   D = LOAD_DATA(matfile)    loads the given .mat file
%
% Returns struct D with fields:
%   D.N        - number of time steps (113)
%   D.Ts       - sample interval [s] (0.036)
%   D.bs       - 3x1 base-station position [x;y;z] (m)
%   D.ant_idx  - 1xN active array index per step (1 or 2)
%   D.ant_or   - 2x2 boresight per array, row k = [az0_k, el0_k] (deg)
%   D.ue_pos   - 3xN TRUE vehicle position (m)  [ERROR ANALYSIS ONLY]
%   D.toa      - 1xN cell, TRUE ranges (m)       [ERROR ANALYSIS ONLY]
%   D.toa_hat  - 1xN cell, MEASURED ranges (m)   [USE THESE]
%   D.aoa      - 1xN cell, TRUE angles 2xM [az;el] (deg) [ERROR ANALYSIS ONLY]
%   D.aoa_hat  - 1xN cell, MEASURED angles 2xm (deg)     [USE THESE]
%   D.z_known  - assumed constant vehicle height (m), 1.5
%
% Pure function: no globals, no plotting.

    if nargin < 1 || isempty(matfile)
        here = fileparts(mfilename('fullpath'));
        matfile = fullfile(here, '..', 'data', 'LNSM_project_data_2026.mat');
    end
    assert(isfile(matfile), 'load_data:missingFile', ...
        'Data file not found: %s', matfile);

    S = load(matfile);

    req = {'ue_pos','bs_pos','antenna_array','antenna_orientation', ...
           'toa','toa_hat','aoa','aoa_hat'};
    for i = 1:numel(req)
        assert(isfield(S, req{i}), 'load_data:missingVar', ...
            'Variable "%s" missing from %s', req{i}, matfile);
    end

    D.N       = size(S.ue_pos, 2);
    D.Ts      = 0.036;                       % from spec (113 steps)
    D.bs      = double(S.bs_pos(:));
    D.ant_idx = double(S.antenna_array(:)'); % 1xN row
    D.ant_or  = double(S.antenna_orientation);
    D.ue_pos  = double(S.ue_pos);
    D.toa     = S.toa;
    D.toa_hat = S.toa_hat;
    D.aoa     = S.aoa;
    D.aoa_hat = S.aoa_hat;
    D.z_known = 1.5;

    % shape checks
    assert(isequal(size(D.bs),     [3 1]), 'load_data:bsShape',  'bs_pos must be 3x1');
    assert(isequal(size(D.ant_or), [2 2]), 'load_data:orShape',  'antenna_orientation must be 2x2');
    assert(numel(D.ant_idx) == D.N,        'load_data:antLen',   'antenna_array length mismatch');
    assert(numel(D.toa)     == D.N && numel(D.toa_hat) == D.N, 'load_data:toaLen', 'toa cell length mismatch');
    assert(numel(D.aoa)     == D.N && numel(D.aoa_hat) == D.N, 'load_data:aoaLen', 'aoa cell length mismatch');
end
