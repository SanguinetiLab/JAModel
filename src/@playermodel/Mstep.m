function [logLik,model] = Mstep(model,U,Y,xhat_plus,Phat_plus,Phatcross_plus,lambda,resp);

T = length(U);
xsz = length(m.x0);

% task parameters
nt = get(m.task,'nt');
Rii = get(m.task,'Rii');
Ri_i = get(m.task,'Ri_i');
ri = get(m.task,'ri');

% Controller model initialization:
if nt==1
    iRii = inv(Rii);
    Hc = -iRii*Ri_i;
    hc = -iRii*ri;
end

% Partner model initialization
P{1} = m.P0;
lambda(1) = m.lambda0;
for t=1:T
   % define partner model parameters
   % Kalman gain (of the partner model):
   K{t} = P{t}*m.C*inv(m.C*P{t}*m.C' + m.SigmaY);
   F{t} = m.A*(eye(length(xsz))-K{t}*m.C);
   G{t} = m.A*K{t}*m.C;
   Sigmaw{t} = m.A*K{t}*m.SigmaY*K{t}'*m.A';
   Sigmav{t} = lambda(t)*iRii;
   
   P{t+1} = F{t}*P{t}*m.A' + m.SigmaX;
   lambda(t+1) = m.a*lambda(t);
end

logLik = 0;
for t=1:T
   B = Hc*Phat_plus{t}*Hc' + Sigmav{t};
   logLik = logLik - 0.5*log(det(B)) ...
       -0.5*(Y(t,:)'-Hc*xhat_plus(:,t) -hc)'*inv(B)*(Y(t,:)'-Hc*xhat_plus(:,t) -hc)';
    
end