function y = fast_log_mvnpdf(X, Mu, Sigma)
%FAST_LOG_MVNPDF Runs fast, no-error-checking version of log(mvnpdf()) for only one data
%point

% X and Mu should be column vectors

% Get size of data.
[d,~] = size(X);

% Set zero mean
X0 = X - Mu;

% Decompose variance
R = chol(Sigma);
xRinv = R' \ X0;
logSqrtDetSigma = sum(log(diag(R)));

% The quadratic form is the inner products of the standardized data
quadform = sum(xRinv.^2, 1);

y = -0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2;
