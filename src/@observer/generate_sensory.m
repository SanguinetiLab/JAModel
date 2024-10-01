function y = generate_sensory(e,u_other,u_self,mode)
switch mode
    case 'sim'

        % y = e.H*u + sqrt(mvnrnd(e.H*u,e.SigmaY,1));
         y = mvnrnd(e.C*u_other + e.D*u_self,e.SigmaY,1);%general ok
        %y = [mvnrnd(-e.D*u_self,e.SigmaY_vis,1) mvnrnd(e.C*u_other + e.D*u_self,e.SigmaY,1) mvnrnd(e.C*u_other,e.SigmaY_vis,1)]; 
       % y = mvnrnd(e.C*u_other,e.SigmaY,1)  + e.D*u_self; %PD cecilia
         
    case 'fit'
        y = e.C*u_other + e.D*u_self;% + sqrt(mvnrnd(e.H*u,e.SigmaY,1));
        
    otherwise
        disp([mode ' is not a recognized modality']);
end
% y = e.H*u + mvnrnd(e.H*u,e.SigmaY,1);