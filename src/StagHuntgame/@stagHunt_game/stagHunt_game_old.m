function g = stagHunt_game(uR,uS,zR,zS,u01,u02,wR,wS)

g.uR = uR;
g.uS = uS;
g.zR = zR;
g.zS = zS;
g.u01 = u01;
g.u02 = u02;

% Cost function weights (Bryson's rule)
g.wR = wR;
g.wS = wS;

% s_i = 1 rabbit strategy
% s_i = 0 stag strategy
s1 = [1 0];
s2 = [1 0];

for s = 1:length(s1)
    % Player 1 strategy 0: Stag
    %           strategy 1: Rabbit
    Rii_1{s} = s1(s)*g.wR + (1-s1(s))*g.wS; 
    Ri_i_1{s} = 0;
    R_i_i_1{s} = (1-s1(s))*g.wS;

    ri_1{s} = - s1(s)*g.wR*(g.uR) - (1-s1(s))*g.wS*(g.uS);
    r_i_1{s} = -(1-s1(s))*g.wS*g.uS;
    z_1{s} = 0.5*s1(s)*g.wR*(g.uR^2 + g.zR) + ...
             0.5*(1-s1(s))*g.wS*(2*g.uS^2 + g.zS);

    % Player 2 strategy 0: Stag
    %           strategy 1: Rabbit
    Rii_2{s} = s2(s)*g.wR + (1-s2(s))*g.wS; 
    Ri_i_2{s} = 0;
    R_i_i_2{s} = (1-s2(s))*g.wS;

    ri_2{s} = - s2(s)*g.wR*(g.uR) - (1-s2(s))*g.wS*(g.uS);
    r_i_2{s} = - (1-s2(s))*g.wS*g.uS;
    z_2{s} = 0.5*s2(s)*g.wR*(g.uR^2 + g.zR) + ...
             0.5*(1-s2(s))*g.wS*(2*g.uS^2 + g.zS);
   
end

g.task1 = task(Rii_1,Ri_i_1,R_i_i_1,ri_1,r_i_1,z_1);
g.task2 = task(Rii_2,Ri_i_2,R_i_i_2,ri_2,r_i_2,z_2);

g = class(g,'stagHunt_game');





