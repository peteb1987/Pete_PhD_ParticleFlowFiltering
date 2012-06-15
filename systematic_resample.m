function [ parents ] = systematic_resample( weights, Np_out )
%SYSTEMATIC_RESAMPLE Systematically resamples a set of weights and returns
%the parent particle for each child.

% weights are assumed to be log-probabilities
weights = exp(weights);

% Number of particles in
Np_in = length(weights);

if nargin < 2
    Np_out = Np_in;
end

% Generate random index array
u = (1/Np_out)*ones(Np_out, 1);
u(1) = rand/Np_out;
u = cumsum(u);

% Generate cumulative weight array
w_sum = cumsum(weights);

% Create array of parent indexes
parents = zeros(Np_out, 1);

% Enumerate offspring
Nchild = zeros(Np_in,1);
Nchild(1) = sum(u < w_sum(1));
parents(1:Nchild(1)) = 1;
cnt = Nchild(1);
for ii = 2:Np_in
    Nchild(ii) = sum((u < w_sum(ii))&(u > w_sum(ii-1)));
    parents(cnt+1:cnt+Nchild(ii)) = ii;
    cnt = cnt + Nchild(ii);
end

end

