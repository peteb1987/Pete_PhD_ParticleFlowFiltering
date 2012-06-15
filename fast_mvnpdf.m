function y = fast_mvnpdf(X, Mu, Sigma)
%MVNPDF Runs fast, no-error-checking version of mvnpdf for only one data
%point

% Get size of data.  Column vectors provisionally interpreted as multiple scalar data.
[n,d] = size(X);

% Set zero mean
X0 = X - Mu;

% Decompose variance
R = chol(Sigma);
xRinv = X0 / R;
logSqrtDetSigma = sum(log(diag(R)));

% The quadratic form is the inner products of the standardized data
quadform = sum(xRinv.^2, 2);

y = exp(-0.5*quadform - logSqrtDetSigma - d*log(2*pi)/2);
