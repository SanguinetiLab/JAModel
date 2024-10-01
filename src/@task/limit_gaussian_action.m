function u = limit_gaussian_action(task,u_i,umin,umax,lambda,Mu,Sigma)
debug = 0;
%% Checking boundaries
H = task.Rii;
% det(H)
f = task.ri + task.Ri_i*u_i;


[u_minc,minval] = quadprog(H,f,[],[],[],[],umin,umax);
[u_maxc,maxval] = quadprog(-H,-f,[],[],[],[],umin,umax);
maxval = -maxval;

if debug
    fplot(@(x) 1/2*x'*H*x+f'*x,[umin umax],'b')
    hold on
    % %fplot(@(x) -1/2*x'*H*x-f'*x,[umin umax],'r')
    fplot(@(x) exp((-1/2*x'*H*x-f'*x)/lambda),[umin umax],'g')
    legend({'Costs','P(action)'})
end

%% Generating action within the desired boundaries
notok = true;
n = length(task.Rii);
mincost = cost(task,u_minc,u_i); % cost for the low cost function and the partner model
p_umin = exp(-mincost/lambda) % max of distribution of actions
maxcost = cost(task,u_maxc,u_i);
p_umax = exp(-maxcost/lambda)

while notok
    u_tent = umin+rand(1,n).*(umax-umin); %mvnrnd(Mu,Sigma,1);%
    p_utent = exp(-cost(task,u_tent,u_i)/lambda);
    rnd = rand(1)*p_umin;
    notok = rnd>p_utent;
    u = u_tent;
end

end