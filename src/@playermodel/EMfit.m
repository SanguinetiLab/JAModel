function model = EMfit(model,U,Y)
% EM fit of player model
loglik_prev = inf;
loglik = 0;
tol = 1e-4;

while abs(loglik-loglik_prev)>tol  
    loglik_prev = loglik;

    % E-step:
   [xhat_plus,Phat_plus,Phatcross_plus,lambda,resp] = Estep(model,U,Y);

    % M-step
   [logLik,model] = ...
       Mstep(model,U,Y,xhat_plus,Phat_plus,Phatcross_plus,lambda,resp);

end