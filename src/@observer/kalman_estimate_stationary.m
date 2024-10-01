function [x_post,P_post,x_prior_next,P_prior_next,K] = kalman_estimate_stationary(e,x_prior,P_prior,y,K_stat)

% Kalman forward single STEP:return current posterior state estimate, and
% nex prior esimate

% Innovation step: combine current prior and current measure to estimate
% current posterior
y_hat = e.C*x_prior;
residual = y - y_hat;
%K = P_prior*e.H'*inv(e.H*P_prior*e.H' + e.SigmaY);
x_post = x_prior + K_stat*residual;
P_post = (eye(e.xsize) - K_stat*e.C)*P_prior;

% Prediction step: use current posterior to estimate next prior
x_prior_next = e.A*x_post;
P_prior_next = e.A*P_post*e.A' + e.SigmaX;




