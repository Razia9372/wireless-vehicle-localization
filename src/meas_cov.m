function R = meas_cov(mode, noise)
%MEAS_COV  Measurement covariance R used by the estimators (NLS & EKF).
%
%   R = meas_cov(mode, noise)
%
% The raw per-quantity noise estimated in Task 1 (sigma_r ~ 0.08 m,
% sigma_az/el ~ 0.5 deg) is the measured-vs-true-LOS-ray scatter. For
% POSITIONING the effective error is larger: the selected min-range ray is
% only an approximate LOS (its range differs from the true geometric
% distance by up to ~0.4 m) and residual multipath leaks in. We therefore
% apply physically-motivated floors so the filter/Gauss-Newton weighting
% and the EKF gating are robust:
%
%   sigma_r  >= 0.20 m,  sigma_az >= 1.0 deg,  sigma_el >= 1.0 deg
%
% Returns R sized for the mode: toa 1x1, aoa 2x2, joint 3x3.
% Pure function.

    sr  = max(noise.sigma_r,  0.20);
    saz = max(noise.sigma_az, 1.00);
    sel = max(noise.sigma_el, 1.00);

    switch lower(mode)
        case 'toa',   R = sr^2;
        case 'aoa',   R = diag([saz^2, sel^2]);
        case 'joint', R = diag([sr^2, saz^2, sel^2]);
        otherwise, error('meas_cov:mode','mode must be toa|aoa|joint');
    end
end
