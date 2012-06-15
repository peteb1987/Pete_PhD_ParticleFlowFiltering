function [ pdf ] = multimodalexp( x, m, P )
%MULTIMODALEXP A multimodal version of the multivariate normal

% x - Nxd - points at which to evaluate pdf
% m - kxd - means
% P - dxdxk - covariances

[N, d] = size(x);
pdf = zeros(N,1);

K = size(m, 1);
assert(size(P,3)==K);

% Loop through data points
for ii = 1:N
    
    logarg = 0;
    
    % Loop through K components
    for kk = 1:K
        
        dist = x(ii,:) - m(kk,:);
        logarg = logarg + log( dist*P(:,:,kk)*dist' );
        
    end
    
    pdf(ii) = exp(-0.5*( exp(logarg/K) ));
    
end


end


