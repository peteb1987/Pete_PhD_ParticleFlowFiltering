function [ y, obs_prb ] = ar1_obs( obs_var, x, y )
%TRACKING_OBS Observation probability and sampling for a direct linear-Gaussian process

if nargin < 3
    y = [];
end

% Sample (if not given)
if isempty(y)
    y = mvnrnd(x, obs_var);
end

% Probability
obs_prb = fast_log_mvnpdf(y, x, obs_var);

end

