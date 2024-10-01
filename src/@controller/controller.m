function c = controller(task,lambda0,a,umin,umax)
% CONTROLLER creates a stochastic controller for a given task
%   Detailed explanation goes here

if nargin==3
    c.lambda0 = lambda0;
    c.a = a;
    c.task = task;
    c.prior = [];
    Rii = get(task, 'Rii');
    if iscell(Rii)
        sz = size(Rii{1},1);
    else
        sz = size(Rii,1);
    end
    c.umin = -inf*ones(sz,1);
    c.umax = +inf*ones(sz,1);

else
    c.lambda0 = lambda0;
    c.a = a;
    c.task = task;
    c.prior = [];
    c.umin = umin;
    c.umax = umax;

end

c = class(c,'controller');

end

