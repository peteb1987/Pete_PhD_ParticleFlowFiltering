function [ pts_array ] = particle_flow_filter( init_pts, times, observs, params, h_trans, h_obs )
%PARTICLE_FLOW_FILTER Run a basic particle flow filter for an AR(1) process

fprintf(1, '\n\n*** Running particle filter. ***\n');

% Initialise constants
Np = length(init_pts);
K = size(times, 2);

% Initialise arrays
pts_array = zeros(Np,K);
last_pts = init_pts;

% Loop through time
for kk = 1:K
    
    fprintf(1, 'Now processing frame %u.\n', kk);
    
    % Create a new particle array
    pts = zeros(Np,1);
    
    % Loop through particles
    for ii = 1:Np
        
        % Propose a new value for the particle
        [pts(ii,1), ~] = feval(h_trans, params.decay, params.proc_var, last_pts(ii,1));
        
    end
    
    %%% Particle flow it to the posterior %%%
    
    % Find Gaussian mean and covariance of the predicted distribution
    m = mean(pts);
    P = var(pts);
    
    H = 1; R = params.obs_var;
    y = observs(1,kk);
    
    % Loop through particles
    for ii = 1:Np
        
        x = pts(ii,1);
        
        dl = params.dl;
        for ll = 0:dl:1
            
            A = -0.5*P*H'*((R+ll*H*P*H')\H);
            b = (1+2*ll*A)*((1+ll*A)*P*H'*(R\y)+A*m);
            x = x + dl*(A*x+b);
            
        end
        
        pts(ii,1) = x;
        
    end
    
    %%% End of particle flow bit %%%
    
    % Store particles and weights
    pts_array(:,kk) = pts;
    
    last_pts = pts;
    
end

end

