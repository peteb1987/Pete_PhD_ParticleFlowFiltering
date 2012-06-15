function ess = ESS( w )
%ESS Calculates effective sample size from array of log-weights

% Convert to linear
w = exp(w);

% Calculate ESS
ess = 1/( sum( w.^2 ) );

end

