function [F,G,Sigma_eta,Sigma_omega] = get_mats(params,e,task,u_other,mode)

for t = 1:T
    y(:,t) = generate_sensory(e,u_other(:,t),mode);
    %[x_post(:,t),P_post{t},x_prior(:,t+1),P_prior{t+1},K{t}] = kalman_estimate(e,x_prior(:,t),P_prior{t},y(:,t));
    [x_post(:,t),P_post{t},x_prior(:,t+1),P_prior{t+1},K{t}] = kalman_estimate(params,x_prior(:,t),P_prior{t},y(:,t));
    
    F{t} = A*(eye(e.xsize) - K{t}*e.H);
    G{t} = A*eye(e.xsize)*e.H;
    Sigma_eta = lambda*inv(task.Rii);
    Sigma_omega{t} = A*K{t}*SigmaY*K{t}'*A';
end
end

