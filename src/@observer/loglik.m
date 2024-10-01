function ll = loglik(e,exp_x,exp_P,exp_crossP,y)

T = length(exp_x);

% initial state estimate 
ll_init_est = -0.5*log(det(e.P0)) - 0.5*trace(inv(e.P0)*(exp_P{1} + (exp_x(:,1) - e.x0)*(exp_x(:,1) - e.x0)'));

% sensory
ll_sensory_t = zeros(e.ysize);
for t = 1:T
    ll_sensory_t = ll_sensory_t + (y(:,t) - e.H*exp_x(:,t))*(y(:,t) - e.H*exp_x(:,t))' + e.H*exp_P{t}*e.H';
end
ll_sens = - (T/2)*log(det(e.SigmaY)) - 0.5*trace(inv(e.SigmaY)*ll_sensory_t);


% estimate
B = zeros(e.xsize); C = zeros(e.xsize); D = zeros(e.xsize);
for t = 2:T
    B = B + (exp_P{t} + exp_x(:,t)*exp_x(:,t)');
    C = C + exp_crossP{t} + exp_x(:,t)*exp_x(:,t-1)';
    D = D + exp_P{t-1} + exp_x(:,t-1)*exp_x(:,t-1)';
end
ll_est = -(T/2)*log(det(e.SigmaX)) - 0.5*trace(inv(e.SigmaX)*(B - C*e.A' - e.A*C' + e.A*D*e.A'));

% complete pseudo - log - likelihood
ll = ll_init_est + ll_sens + ll_est;

