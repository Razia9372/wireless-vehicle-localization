function noise = load_noise(root)
%LOAD_NOISE  Load the Task-1 measurement-noise model used by the estimators.
%
%   noise = load_noise(root)   root = project root folder
%
% The noise model (sigma_r, sigma_az, sigma_el and the R matrices) is
% estimated from the measured-vs-true LOS residuals by
% scripts/task1_data_analysis.m and saved to results/noise_model.mat.
% Fails with a clear message if Task 1 has not been run yet.

    f = fullfile(root, 'results', 'noise_model.mat');
    if ~isfile(f)
        error('load_noise:missing', ...
            ['results/noise_model.mat not found.\n' ...
             'Run scripts/task1_data_analysis.m (or scripts/run_all.m) first - ' ...
             'it estimates the measurement noise the estimators need.']);
    end
    L = load(f);
    noise = L.noise;
end
