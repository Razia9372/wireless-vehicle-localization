%RUN_ALL  Regenerate the entire LNSM project end-to-end.
%
% Runs Tasks 1-4 in order and (optionally) the tracking animation, saving
% every figure to results/figures/ and the metrics table to results/.
% Each task runs in an isolated workspace so its `clear` does not affect the
% others.
%
% Usage:  run('scripts/run_all.m')      from the project root, or just run
%         this file from the scripts/ folder.

function run_all()
    here = fileparts(mfilename('fullpath'));
    root = fullfile(here,'..');
    addpath(fullfile(root,'src'));
    use_light_figures();
    t0 = tic;

    fprintf('\n############## TASK 1: data analysis ##############\n');
    run_isolated(fullfile(here,'task1_data_analysis.m'));
    fprintf('\n############## TASK 2: NLS vs EKF ##############\n');
    run_isolated(fullfile(here,'task2_nls_vs_ekf.m'));
    fprintf('\n############## TASK 3: motion models ##############\n');
    run_isolated(fullfile(here,'task3_motion_models.m'));
    fprintf('\n############## TASK 4: evaluation ##############\n');
    run_isolated(fullfile(here,'task4_evaluation.m'));

    fprintf('\n############## (optional) tracking animation ##############\n');
    try
        make_tracking_video(root);
    catch ME
        fprintf('Animation skipped (%s)\n', ME.message);
    end

    fprintf('\nAll done in %.1f s. Figures in results/figures, video in results/video.\n', toc(t0));
end

function run_isolated(scriptfile)
    run(scriptfile);              % isolated: this function's workspace
end

function make_tracking_video(root)
    addpath(fullfile(root,'src'));
    D = load_data();
    noise = load_noise(root);
    est = ekf_track(D,'cv','joint',noise);
    out = make_animation(D, est, 'EKF-cv joint', ...
        fullfile(root,'results','video','tracking_ekf_cv_joint'), struct('fps',15));
    fprintf('Wrote %s\n', out);
end
