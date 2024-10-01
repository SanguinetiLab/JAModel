function u = rndgibbs(task,u_minc,u_i,umin,umax,lambda)
% RNDGIBBS generate a random action in Gibbs case.
%   umin and umax are the minumum and maximum values for the action.

notok = true;
n = length(task.Rii);
mincost = cost(task,u_minc,u_i); % cost for the low cost function and the partner model
p_umax = exp(-mincost/lambda); % max of distribution of actions 

while notok
    u_tent = umin+rand(1,n).*(umax-umin);
    p_utent = exp(-cost(task,u_tent,u_i)/lambda);
    rnd = rand(1)*p_umax;
    notok = rnd>p_utent;
    u = u_tent;
end

end

