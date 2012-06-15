function [ x, trans_prb ] = ar1_trans( decay, proc_var, last_x, x )
%TRACKING_TRANS Transition probability and sampling for an AR(1) process

if nargin < 4
    x = [];
end

% Prediction mean
x_mn = decay*last_x;

% Sample (if not given)
if isempty(x)
    x = mvnrnd(x_mn, proc_var);
end

% Probability
trans_prb = fast_log_mvnpdf(x, x_mn, proc_var);

end

