function [xhat_plus,Phat_plus,Phatcross_plus,lambda,resp] = Estep(m,U,Y);

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



% Partner model + Controller is defined by:
% xi(t+1) = F(t) xi(t) + G(t)*u-i(t) + w(t)
% y(t)=ui(t) = Hc*xi(t) + hc + v(t)

% E-step requires P+C model in predictor form:
W = cell(T,1);
xhat_plus=zeros(1,T);
xhat_minus=zeros(1,T);
Phat_minus = cell(T,1);
Phat_plus = cell(T,1);
Phatcross_plus = cell(T,1);
Phat_minus{1}=m.P0;
xhat_minus(:,1)=m.x0;
resp = zeros(T,nt);

% forward step
for t=1:T
   % Measurement update (correction step)
   W{t} = Phat_minus{t}*Hc*inv(Hc*Phat_minus{t}*Hc' + Sigmav{t});
   xhat_plus(:,t) = xhat_minus(:,t)+W{t}*(Y(t,:)'-Hc*xhat_minus(:,t)-hc);
   Phat_plus{t} = (eye(xsz)-W{t}*Hc)*Phat_minus{t};

   % Time update (prediction step)
   xhat_minus(:,t+1) = F{t}*xhat_plus(:,t) + G{t}*U(t,:)'; 
   Phat_minus{t+1} = F{t}*Phat_plus{t}*F{t}' + Sigmaw{t};   
end

% backward step
xxhat_plus(:,T)=xhat_plus(:,T);
PPhat_plus{T} = Phat_plus{T};
for t=T:-1:2
  J{t-1} = Phat_plus{t-1}*F{t-1}*inv(Phat_plus{t-1});
  xxhat_plus(:,t-1) = xhat_plus(:,t-1) + ...
                      J{t-1}*(xxhat_plus(:,t)-F{t-1}*xhat_plus(:,t-1));
  PPhat_plus{t-1} = Phat_plus{t-1} + ...
                      J{t-1}*(PPhat_plus{t}-Phat_minus{t})*J{t-1}';
end

Vhatcross_plus{T} = (eye(xsz)-W{T}*Hc)*F{T}*Phat_plus{T-1};
Phatcross_plus{T}=Vhatcross_plus{T}+xxhat_plus(:,T)*xxhat_plus(:,T)'; 
for t=T:-1:3
  Vhatcross_plus{t-1}=...'
      Phat_plus{t-1}*J{t-2}'+...
      J{t-1}*(Vhatcross_plus{t}-F{t-1}*Phat_plus{t-1})*J{t-2}';
  Phatcross_plus{t-1}=Vhatcross_plus{t-1}+...
      xxhat_plus(:,t-1)*xxhat_plus(:,t-1)';
end
xhat_plus = xxhat_plus;
Phat_plus = PPhat_plus;
