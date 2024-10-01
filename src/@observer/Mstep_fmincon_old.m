function obs = Mstep_fmincon(e,exp_x,exp_P,exp_crossP,y)

T = length(exp_x);

B = zeros(e.xsize); C = zeros(e.xsize); D = e.xsize; E = e.ysize;

for t = 2:T
    B = B + exp_P{t-1} + exp_x(:,t-1)*exp_x(:,t-1)';
    C = C + exp_crossP{t-1} + exp_x(:,t)*exp_x(:,t-1)';
    D = D + exp_P{t}  + exp_x(:,t)*exp_x(:,t)';
    E = E + e.H*((y(:,t) - exp_x(:,t))*(y(:,t) - exp_x(:,t))' + exp_P{t})*e.H';
end

params0 = [e.A zeros(e.xsize,2*e.xsize + e.ysize + 1);...
    zeros(e.xsize) e.SigmaX zeros(e.xsize,e.ysize + 1 + e.xsize); ...
    zeros(e.ysize,2*e.xsize) e.SigmaY zeros(e.ysize, 1 + e.xsize);...
    zeros(e.xsize,2*e.xsize + e.ysize) e.x0 zeros(e.xsize);...
    zeros(e.xsize,2*e.xsize + e.ysize + 1) e.P0];

params_low = [eye(length(e.A)).*eps zeros(e.xsize,2*e.xsize + e.ysize + 1);...
    zeros(e.xsize) eye(length(e.SigmaX))*eps zeros(e.xsize,e.ysize + 1 + e.xsize); ...
    zeros(e.ysize,2*e.xsize) eye(length(e.SigmaY))*eps zeros(e.ysize, 1 + e.xsize);...
    zeros(e.xsize,2*e.xsize + e.ysize) zeros(size(e.x0)) zeros(e.xsize);...
    zeros(e.xsize,2*e.xsize + e.ysize + 1) eye(length(e.P0)).*eps];

params_up= [eye(length(e.A)).*1 zeros(e.xsize,2*e.xsize + e.ysize + 1);...
    zeros(e.xsize) eye(length(e.SigmaX))*100 zeros(e.xsize,e.ysize + 1 + e.xsize); ...
    zeros(e.ysize,2*e.xsize) eye(length(e.SigmaY))*100 zeros(e.ysize, 1 + e.xsize);...
    zeros(e.xsize,2*e.xsize + e.ysize) ones(size(e.x0)).*100 zeros(e.xsize);...
    zeros(e.xsize,2*e.xsize + e.ysize + 1) eye(length(e.P0)).*100];

Apos = [1 e.xsize 1 e.xsize];
SigX_pos = [Apos(2)+1  Apos(2)+e.xsize Apos(4)+1 Apos(4)+e.xsize] ;
SigY_pos = [SigX_pos(2)+1  SigX_pos(2)+e.ysize SigX_pos(4)+1 SigX_pos(4)+e.ysize];
x0_pos = [SigY_pos(2)+1  SigY_pos(2)+e.xsize SigY_pos(4)+1 SigY_pos(4)+1];
P0_pos = [x0_pos(2)+1  x0_pos(2)+e.xsize x0_pos(4)+1 x0_pos(4)+e.xsize];

fun = @(params) -(-0.5*log(det(params(P0_pos(1):P0_pos(2),P0_pos(3):P0_pos(4)))) ...
    - 0.5*trace(inv(params(P0_pos(1):P0_pos(2),P0_pos(3):P0_pos(4)))*(exp_P{1} + (exp_x(:,1) - params(x0_pos(1):x0_pos(2),x0_pos(3):x0_pos(4)))*(exp_x(:,1) - params(x0_pos(1):x0_pos(2),x0_pos(3):x0_pos(4)))'))...
    - (T/2)*log(det(params(SigX_pos(1):SigX_pos(2),SigX_pos(3):SigX_pos(4))))...
    - 0.5*trace(inv(params(SigX_pos(1):SigX_pos(2),SigX_pos(3):SigX_pos(4)))*(B - C*params(Apos(1):Apos(2),Apos(3):Apos(4))' - params(Apos(1):Apos(2),Apos(3):Apos(4))*C' + params(Apos(1):Apos(2),Apos(3):Apos(4))*D*params(Apos(1):Apos(2),Apos(3):Apos(4))'))...
    - (T/2)*log(det(params(SigY_pos(1):SigY_pos(2),SigY_pos(3):SigY_pos(4))))...
    - 0.5*trace(inv(params(SigY_pos(1):SigY_pos(2),SigY_pos(3):SigY_pos(4)))*E));


[params_new,fval,exitflag,output] = fmincon(fun, params0,[],[],[],[],params_low,params_up);

A_new = params_new(Apos(1):Apos(2),Apos(3):Apos(4));
SigmaX_new = params_new(SigX_pos(1):SigX_pos(2),SigX_pos(3):SigX_pos(4));
SigmaY_new = params_new(SigY_pos(1):SigY_pos(2),SigY_pos(3):SigY_pos(4));
x0_new = params_new(x0_pos(1):x0_pos(2),x0_pos(3):x0_pos(4));
P0_new = params_new(P0_pos(1):P0_pos(2),P0_pos(3):P0_pos(4));

obs = observer(A_new,SigmaX_new,SigmaY_new,eye(e.ysize),x0_new,P0_new);

end

