function u_minc = maxgibbs(task,u_i,umin,umax) 
% MAXGIBBS calculate the min value of the cost function
% written as a quadratic form - function name is misleading since it
% is not just for gibbs distribution

H = task.Rii;
f = task.ri + task.Ri_i*u_i;

%fplot(@(x) 1/2*x'*H*x+f'*x,[umin umax],'b')

[u_minc,minval] = quadprog(H,f,[],[],[],[],umin,umax);

% [u_maxc,maxval] = quadprog(-H,-f,[],[],[],[],umin,umax);
% maxval = -maxval;
% hold on
% fplot(@(x) -1/2*x'*H*x-f'*x,[umin umax],'r')
end

