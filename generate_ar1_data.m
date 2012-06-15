function [ t, x, y ] = generate_ar1_data( params )
%GENERATE_AR1_DATA Generate AR(1) data using exact parameters

% Create arrays
t = cumsum(ones(1,params.K))*params.dt;
x = zeros(1,params.K);
y = zeros(1,params.K);

% Initialise state
last_x = params.x0;

% Sampling recursion
for kk = 1:params.K
    
    % Sample state
    [x(kk), ~] = ar1_trans(params.decay, params.proc_var, last_x);
    
    % Sample observation
    [y(kk), ~] = ar1_obs(params.obs_var, x(kk));
    
    last_x = x(kk);
    
end

end

