% Set parameters

global params;

%Model
params.dt = 1;              % Time step size
params.K = 1000;             % Number of time steps
params.proc_var = 1;        % Process variance
params.obs_var = 10;         % Observation variance
params.decay = 0.95;        % Decay coefficient
params.x0 = 0;              % Starting value

% Priors (hyperparameters) : Normal for starting state. Normal-Gamma for decay and process precision
hyper_params.start_prior_mn = 1;          % Starting state prior: mean
hyper_params.start_prior_var = 1;         % Starting state prior: variance
hyper_params.decay_prior_mn = 0;          % Decay prior: mean
hyper_params.decay_prior_df = 1;          % Decay prior: degrees of freedom (for normal-gamma prior)
hyper_params.procprec_prior_shape = 2;    % Process precision prior: shape
hyper_params.procprec_prior_scale = 0.5;  % Process precision prior: scale

% Algorithm
params.Nx = 100;
params.dl = 0.1;

