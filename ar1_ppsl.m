function [ x, ppsl_prb ] = ar1_ppsl( decay, proc_var, obs_var, last_x, y, x)
%TRACKING_PPSL Propose a new state for an AR(1) process

if nargin < 6
    x = [];
end

% Kalman filter gives optimal prediction
[m, P] = kf_predict(last_x, 0, decay, proc_var);
[m, P] = kf_update(m, P, y, 1, obs_var);

% Sample (if not given)
if isempty(x)
    x = mvnrnd(m, P);
end

% Probability
ppsl_prb = fast_log_mvnpdf(x, m, P);

end

