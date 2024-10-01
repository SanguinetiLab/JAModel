function [x_post,P_post,x_prior_next,P_prior_next,K] = kalman_estimate(e,x_prior,P_prior,y,u)

% Kalman forward single STEP:return current posterior state estimate, and
% nex prior esimate


% Innovation step: combine current prior and current measure to estimate
% current posterior
if isempty(u)
   y_hat = e.C*x_prior;
else    
   y_hat = e.C*x_prior + e.D*u;
end

residual = y - y_hat;

K  = P_prior*e.C'*inv(e.C*P_prior*e.C' + e.SigmaY);

x_post = x_prior + K*residual;
P_post = (eye(e.xsize) - K*e.C)*P_prior;

% Prediction step: use current posterior to estimate next prior
if isempty(u)
    x_prior_next = e.A*x_post;
else
    x_prior_next = e.A*x_post;% + e.B*u;
end
P_prior_next = e.A*P_post*e.A' + e.SigmaX;




