function e = Mstep(e,exp_x,exp_P,exp_crossP,y,u_other,equation)
% function e = Mstep(e,exp_x,exp_P,y,equation)
T = length(exp_x);

B = 0; C = 0; D = 0;
for t = 2:T
    B = B + exp_P{t-1} + exp_x(:,t-1)*exp_x(:,t-1)';
%     C = C + exp_crossP{t-1} + exp_x(:,t)*exp_x(:,t-1)';
    D = D + exp_P{t-1}  + exp_x(:,t)*exp_x(:,t)';
end

xsize = get(e,'xsize');
ysize = get(e,'ysize');

ssy = zeros(ysize);
switch equation
    case 'Sanguineti_DeVicariis'
        a = trace(C+C')/trace(B); % from Cecilia computation - considering scalar parameters
        e.A = 0.8;%eye(xsize).*a;%
        
       
        sigmax = 0.001^2;%trace(A^2.*B - A.*(C + C') + D)/((T-1));%
        e.SigmaX = eye(xsize).*sigmax;
       
        for t = 1:T
            ssy = ssy + ...
                e.H*((y(:,t) - exp_x(:,t))*(y(:,t) - exp_x(:,t))'+ exp_P{t})*e.H';
        end
        
        sigmay = 0.001^2;%trace(ssy)/T;%
        e.SigmaY = eye(ysize).*sigmay;
        e.P0 = trace(exp_P{1} + exp_x(:,1)*exp_x(:,1)').*eye(xsize); %eps*eye(xsize);%
    case 'Shumway_Stoffer'
        %e.A =C*B';%0.8;%
        
        %e.SigmaX = (1/(T-1))*(D - C*inv(B)*C');%

        for t = 1:T
            ssy = ssy + (y(:,t) - e.H*exp_x(:,t))*(y(:,t) - e.H*exp_x(:,t))' + e.H*exp_P{t}*e.H';
        end
        %e.SigmaY =  0.001^2;%ssy/T;%
        e.P0 = exp_P{1} - exp_x(:,1)*exp_x(:,1)' ;
        
    case 'Hinton'
        a1 = zeros(size(exp_crossP{1})); a2 = a1;
        for t = 2:T
            a1 = a1+exp_crossP{t};
            a2 = a2+exp_P{t};
        end
       %e. A = a1*inv(a2);
        
        ssy = zeros(size(exp_crossP{1}));
        for t = 2:T
%             ssy = ssy + (y(:,t)*y(:,t)' - e.H*exp_x(:,t)*y(:,t)');
            ssy = ssy + e.H*(u_other(:,t) - exp_x(:,t)*u_other(:,t)'*e.H);
        end
        e.SigmaY = (1/T)*ssy;
        
        ssx1 = zeros(size(exp_crossP{1})); ssx2 = ssx1;
        for t = 2:T
           ssx1 = ssx1 + exp_P{t};
           ssx2 = ssx2 + exp_crossP{t} ;
        end
            
        %e.SigmaX = (1/(T-1))*(ssx1 - A*ssx2);
        
        %e.P0 = exp_P{1} - exp_x(:,1)*exp_x(:,1)' ;
end

% e.x0 = exp_x(:,1);%zeros(size(exp_x(:,1)));%


%e.H = eye(xsize);


%e = observer(A,SigmaX,SigmaY,H,x0,P0);

