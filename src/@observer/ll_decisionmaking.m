function logL = ll_decisionmaking(params,e,t,u,u_other,mode)

T = length(exp_x);

[F,G,Sigma_eta,Sigma_omega] = get_mats(params,e,t,T,u_other,mode);

LDS.A = F;              % ndirs x ndirs
LDS.fixA = 0;
LDS.B = G;                             % ndirs x ndirs (but 5 free params)  
LDS.fixB = 0;
LDS.C = eye(e.xsize);                % 1 x ndirs
LDS.fixC = 1;
% LDS.D = 0;%D';                            % ndirs x 1     (but 3 free params)
% LDS.fixD = 0;
LDS.Q = Sigma_omega;           % ndirs x ndirs (1 free param)
LDS.fixQ = 0;
LDS.R = Sigma_eta; %*eye(ndirs,ndirs);           % ndirs x ndirs (1 free param)
LDS.fixR = 1;
LDS.x0 = e.x0;  % ndirs x 1 
LDS.fixx0 = 0;
LDS.V0 = e.P0;    % ndirs x ndirs
LDS.fixV0 = 0;
        
Y = u;        % output: 1 x T
ny = size(Y,1);
       
        
DIR = direction';         % direction: used to determine direction of each individual trial
                          
U = zeros(ndirs,T);        % input (this is actually feedback): ndirs x T
for i = 1:T
  U(DIR(1,i),i) = Y(1,i);   
end
nu = size(U,1);
        
% E-step (initial estimate of state)
[Lik,X,V,V1,SUMS,Yest] = SmoothLDS(LDS,Y,U);
        
logL = -Lik;
end

