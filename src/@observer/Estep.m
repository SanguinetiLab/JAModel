function [x_prior,P_prior,x_post,P_post,exp_x,exp_P,exp_crossP,y] = Estep(e,u_other)
T = length(u_other);
x_prior(:,1) = e.x0;
P_prior{1} = e.P0;

% Forward Iteration
for t = 1:T
    y(:,t) = generate_sensory(e,u_other(:,t),'fit');
    
    [x_post(:,t),P_post{t},x_prior(:,t+1),P_prior{t+1},K{t}] = kalman_estimate(e,x_prior(:,t),P_prior{t},y(:,t));
end

exp_x = zeros(size(y));
exp_P = cell(T,1);
exp_crossP = cell(T,1);
J = cell(T,1);

exp_x(:,T) = x_post(:,T);
exp_P{T} =  P_post{T};

% Backward Iteration
for t = T-1:-1:1 
    J{t} = P_post{t} *e.A'*inv(P_post{t+1}); % 
    exp_x(:,t) = x_post(:,t) + J{t}*(exp_x(:,t+1) - x_post(:,t+1));
    exp_P{t}  = P_post{t}  + J{t}*(exp_P{t+1}  - P_post{t+1})*J{t}';
end
exp_crossP{T}  = (eye(e.xsize)-K{T}*e.H)*e.A*exp_P{T};

% Cross
for t = T-1:-1:1
   % exp_crossP(:,:,t) = P_post(:,:,t)*J(:,:,t-1)' + J(:,:,t)*(exp_crossP(:,:,t) - e.A*P_post(:,:,t))*J(:,:,t-1)';
exp_crossP{t} = exp_P{t+1}*J{t};
end
