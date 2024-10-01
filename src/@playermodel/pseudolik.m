function [logQ,uhat_minus,Phat_plus,pi_posterior] = pseudolik(m,U,Y)
% Y is player's action --> the output
% U is the partner's action --> the input

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
    Hc_hat = -iRii*Ri_i;
    hc_hat = -iRii*ri';
else
    for n=1:nt
        iRii{n} = inv(Rii{n});
        Hc_hat{n} = -iRii{n}*Ri_i{n};
        hc_hat{n} = -iRii{n}*ri{n}';
    end
end

% Partner model initialization
Pm{1} = m.P0; % P^{-}

% Define partner model parameters 
% (calculated once as they do not depend on data)
for t=1:T 
    % State equation of Partner model + Controller (P+C):
    % xi(t+1) = F(t) xi(t) + G(t)*u-i(t) + w(t)
    
    % Kalman gain (of partner model):
    K{t} = Pm{t}*m.C'*inv(m.C*Pm{t}*m.C' + m.SigmaY); 
    
    F{t} = m.A*(eye(xsz)-K{t}*m.C);
    G{t} = m.A*K{t}*m.C; % ok
    Sigmaw{t} = m.A*K{t}*m.SigmaY*K{t}'*m.A'; 

    % Update equation for partner model:
    Pm{t+1} = F{t}*Pm{t}*m.A' + m.SigmaX; 
end

% Partner model + Action controller is defined by:
% xi(t+1) = F(t) xi(t) + G(t)*u-i(t) + w(t)
% ui(t) = Hc*xi(t) + hc + v(t)
%
% where: 
% partner action u-i is system input
% player output ui is system output

% E-step requires P+C model in predictor form:
% 
% Kalman gain: 
% W(t) = Phat_minus(t)*Hc'*inv(Hc*Phat_minus(t)*Hc' + Sigmav(t));
%
% Prediction of output:
% uihat_minus(t) = Hc*xhat_minus(t)+hc
% Qhat_minus(t) = Hc*Phat_minus(t)*Hc'  
%
% Correction:
% xhat_plus(t) =  xhat_minus(t) + W(t) [ui(t)-uihat_minus(t)]
% Phat_plus(t) =  (I-W(t)*Hc)*Phat_minus(t);
%
% Update:
% xhat_minus(t+1) = F(t)*xhat_plus(t) + G(t)*u-i(t)
% Phat_minus(t+1) = F(t)*Phat_plus(t)*F(t)' + Sigmaw(t)

W = cell(T,1);

xhat_plus=zeros(xsz,T);
xhat_minus=zeros(xsz,T);

Phat_minus = cell(T,nt);
Phat_plus = cell(T,1);

uhat_minus = cell(T,nt);%zeros(xsz,T);
Qhat_minus = cell(T,nt);

% Initialization
for n = 1:nt
    Phat_minus{1,n}=m.P0;
end
xhat_minus(:,1)=m.x0;
lambda(1) = m.lambda0;

% forward step
for t=1:T
    if nt == 1

        % Action covariance
        Sigmav{t} = get_covariance(m.task,lambda(t)); % this is lambda*inv(Rii)

        % Predicted output (action)
        uhat_minus{t} = Hc_hat*xhat_minus(:,t)+hc_hat;
        Qhat_minus{t} = Hc_hat*Phat_minus{t}*Hc_hat' + Sigmav{t}; %THIS WAS ERROR!

        % Kalman gain:
        W{t} = Phat_minus{t}*Hc_hat*inv(Qhat_minus{t});%  %YEP consisent with VS notation B = Phat

        % Correction step        
        xhat_plus(:,t) = xhat_minus(:,t)+W{t}*(Y(:,t)'-uhat_minus{t}); % ok, double check but I agree that Y is the output, player action
       
        Phat_plus{t} = (eye(xsz)-W{t}*Hc_hat)*Phat_minus{t}; % OK
        
        pi_posterior(t,:) = 1;
    else
        % Action parameters:
        Mu{t} = get_center(m.task,xhat_minus(:,t)); 
        Sigmav{t} = get_covariance(m.task,lambda(t));
        priors{t} = get_priors(m.task,xhat_minus(:,t),Phat_minus{t},lambda(t));
        
        gm = gmdistribution(Mu{t}',Sigmav{t},priors{t});

        % Need to compute posterior probability of strategy given partner
        % model and observation
        [idx,nlogL,action_prob(t,:)] = cluster(gm,Y(:,t)');

        % for n = 1:nt
        %     pui_sixi(t,n) = mvnpdf(Y(:,t)',Mu{t}(n),Sigmav{t}(:,:,n))
        % end
        pi_posterior(t,:) = action_prob(t,:);
        % pi_posterior is actually equal to action_prob
        
        xhat_plus(:,t) = xhat_minus(:,t);
        Phat_plus{t} = zeros(xsz);  
        
        for n = 1:nt

            % Predicted output (action) - per strategy
            uhat_minus{t,n} = Hc_hat{n}*xhat_minus(:,t)  + hc_hat{n}; %ok
            Qhat_minus{t,n} = Hc_hat{n}*Phat_minus{t}*Hc_hat{n}' + Sigmav{t}(:,:,n); %THIS WAS ERROR!
            
            % Kalman gain - per strategy
            Wk{n} = Phat_minus{t}*Hc_hat{n}*inv(Qhat_minus{t,n});% +Sigmav{t}(:,:,n)); %ok

            % Correction step - weighted by strategy
            xhat_plus(:,t) = xhat_plus(:,t) + pi_posterior(t,n)*Wk{n}*(Y(:,t) - uhat_minus{t,n}); 
            Phat_plus{t} = Phat_plus{t} + pi_posterior(t,n)*(eye(xsz) - Wk{n}*Hc_hat{n})*Phat_minus{t}; %ok
        end
    end

    % Time update (prediction step)
    xhat_minus(:,t+1) = F{t}*xhat_plus(:,t) + G{t}*U(:,t);
    Phat_minus{t+1} = F{t}*Phat_plus{t}*F{t}' + Sigmaw{t};
    
    % Lambda update
    lambda(t+1) = m.a*lambda(t);

end

% backward step
% xxhat_plus(:,T)=xhat_plus(:,T);
% PPhat_plus{T} = Phat_plus{T};
% for t=T:-1:2
%   J{t-1} = Phat_plus{t-1}*F{t-1}*inv(Phat_plus{t-1});
%   xxhat_plus(:,t-1) = xhat_plus(:,t-1) + ...
%                       J{t-1}*(xxhat_plus(:,t)-F{t-1}*xhat_plus(:,t-1));
%   PPhat_plus{t-1} = Phat_plus{t-1} + ...
%                       J{t-1}*(PPhat_plus{t}-Phat_minus{t})*J{t-1}';
% end
%
% Vhatcross_plus{T} = (eye(xsz)-W{T}*Hc)*F{T}*Phat_plus{T-1};
% Phatcross_plus{T}=Vhatcross_plus{T}+xxhat_plus(:,T)*xxhat_plus(:,T)';
% for t=T:-1:3
%   Vhatcross_plus{t-1}=...'
%       Phat_plus{t-1}*J{t-2}'+...
%       J{t-1}*(Vhatcross_plus{t}-F{t-1}*Phat_plus{t-1})*J{t-2}';
%   Phatcross_plus{t-1}=Vhatcross_plus{t-1}+...
%       xxhat_plus(:,t-1)*xxhat_plus(:,t-1)';
% end
% xhat_plus = xxhat_plus;
% Phat_plus = PPhat_plus;

% -E{logLik} in innovation form...
% % From Shumway & Stoffer (1982) taken from Gupta & Mehra (1974)
% logQ = 0;
% for t=1:T
%     logQ = logQ + log(det(Pyhat_minus{t}))+...
%         (Y(t,:)'-yhat_minus(:,t))'*...
%         inv(Pyhat_minus{t})*...
%         (Y(t,:)'-yhat_minus(:,t));
% end
% logQ = 0.5*logQ;
%
% % General formulation:
% logQ = 0;
% for t=1:T
%     logQ = logQ + log(det(Pyhat_minus{t}))+...
%         (Y(t,:)'-yhat_minus(:,t))'*...
%         inv(Pyhat_minus{t})*...
%         (Y(t,:)'-yhat_minus(:,t));
% end
% d = length(yhat_minus(:,1));
% logQ = 0.5*logQ + T*(d/2)*log(2*pi);

% General (multi-strategy) formulation:
logQ = 0;
if nt>1
    for t=1:T
        Qtn = 0;
        for n = 1:nt
            Qtn = Qtn + pi_posterior(t,n)*(1./sqrt(det(Qhat_minus{t,n})))*...
                exp(-0.5*((Y(:,t)-uhat_minus{t,n})'*inv(Qhat_minus{t,n})*...
                          (Y(:,t)-uhat_minus{t,n})));
            sd(n,t) = pi_posterior(t,n)*(1./sqrt(det(Qhat_minus{t,n})));
            du(n,t) = exp(-0.5*((Y(:,t)-uhat_minus{t,n})'*inv(Qhat_minus{t,n})*...
                          (Y(:,t)-uhat_minus{t,n})));
        end
        logQ = logQ + log(Qtn+eps);
    end
    d = length(uhat_minus{1,1});    
    logQ = -logQ;
    logQ = 0.5*logQ + T*(d/2)*log(2*pi);
else
    
    for t=1:T
        logQ = logQ + log(det(Qhat_minus{t}))+...
            (Y(:,t)-uhat_minus{t})'*inv(Qhat_minus{t})*...
            (Y(:,t)-uhat_minus{t});
    end
    d = length(uhat_minus{1});
    logQ = 0.5*logQ + T*(d/2)*log(2*pi);
    
end

