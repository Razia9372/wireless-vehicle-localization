function mm = motion_models(model, Ts, q)
%MOTION_MODELS  Linear state-transition F and process-noise Q (Task 3).
%
%   mm = motion_models(model, Ts, q)
%
% Horizontal-plane kinematic models (z fixed). Returns struct mm:
%   .name    - model id
%   .dim     - state dimension
%   .F       - dim x dim state-transition matrix
%   .Q       - dim x dim process-noise covariance
%   .posIdx  - indices of [x; y] within the state vector
%   .velIdx  - indices of [vx; vy] (empty for CP)
%   .lift    - @(Hp) -> full meas Jacobian: places the 2 position columns
%              of Hp (rows x 2) at posIdx, zeros elsewhere
%   .x0extra - template extra (velocity/accel) entries for initialisation
%
% Models (per-axis blocks, then interleaved x-block then y-block):
%   'cp' constant position / random walk : state [x; y]
%   'cv' constant velocity               : state [x; vx; y; vy]
%   'ca' constant acceleration           : state [x; vx; ax; y; vy; ay]
%
% q is the process-noise intensity (CP: position var; CV: white-noise
% acceleration PSD; CA: white-noise jerk PSD).
%
% Pure function.

    switch lower(model)
        case 'cp'
            F = eye(2);
            Q = q * eye(2);
            mm.dim = 2; mm.posIdx = [1 2]; mm.velIdx = [];

        case 'cv'
            Fb = [1 Ts; 0 1];
            Qb = q * [Ts^3/3, Ts^2/2; Ts^2/2, Ts];
            F  = blkdiag(Fb, Fb);
            Q  = blkdiag(Qb, Qb);
            mm.dim = 4; mm.posIdx = [1 3]; mm.velIdx = [2 4];

        case 'ca'
            Fb = [1 Ts Ts^2/2; 0 1 Ts; 0 0 1];
            Qb = q * [Ts^5/20, Ts^4/8, Ts^3/6;
                      Ts^4/8,  Ts^3/3, Ts^2/2;
                      Ts^3/6,  Ts^2/2, Ts];
            F  = blkdiag(Fb, Fb);
            Q  = blkdiag(Qb, Qb);
            mm.dim = 6; mm.posIdx = [1 4]; mm.velIdx = [2 5];

        otherwise
            error('motion_models:model', 'model must be cp|cv|ca, got "%s"', model);
    end

    mm.name = lower(model);
    mm.F = F; mm.Q = Q;
    pidx = mm.posIdx; dim = mm.dim;
    mm.lift = @(Hp) lift_jac(Hp, pidx, dim);
end

function Hfull = lift_jac(Hp, posIdx, dim)
    Hfull = zeros(size(Hp,1), dim);
    Hfull(:, posIdx) = Hp;
end
