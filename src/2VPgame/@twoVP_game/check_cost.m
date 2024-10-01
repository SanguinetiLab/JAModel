function [Jiok,J1,J2]=check_cost(g,ui,u_i,n)
% COST calculates the cost of a task (part of a quadratic game) for given actions U1,U2
phs = phspline(g.times);

g.giVP{1} = gmtx(phs,g.tcE);
g.giVP{2} = gmtx(phs,g.tcL);
g.giVP{3} = gmtx(phs,g.tcM);

g.GE = g.giVP{1}'*g.giVP{1};
g.GL = g.giVP{2}'*g.giVP{2};
g.GM = g.giVP{3}'*g.giVP{3};

QVP{1} = g.Sx'*g.GE*g.Sx + g.Sy'*g.GE*g.Sy;
QVP{2} = g.Sx'*g.GL*g.Sx + g.Sy'*g.GL*g.Sy;
QVP{3} = g.Sx'*g.GM*g.Sx + g.Sy'*g.GM*g.Sy;

QVP0{1} = g.Sx0'*g.GE*g.Sx + g.Sy0'*g.GE*g.Sy;
QVP0{2} = g.Sx0'*g.GL*g.Sx + g.Sy0'*g.GL*g.Sy;
QVP0{3} = g.Sx0'*g.GM*g.Sx + g.Sy0'*g.GM*g.Sy;

QVP00{1} = g.Sx0'*g.GE*g.Sx0 + g.Sy0'*g.GE*g.Sy0;
QVP00{2} = g.Sx0'*g.GL*g.Sx0 + g.Sy0'*g.GL*g.Sy0;
QVP00{3} = g.Sx0'*g.GM*g.Sx0 + g.Sy0'*g.GM*g.Sy0;

Q = g.Sx'*g.G*g.Sx + g.Sy'*g.G*g.Sy;
Q0 = g.Sx0'*g.G*g.Sx + g.Sy0'*g.G*g.Sy;
Q00 = g.Sx0'*g.G*g.Sx0 + g.Sy0'*g.G*g.Sy0;

Qd = g.Sx'*g.Gd*g.Sx + g.Sy'*g.Gd*g.Sy;
Q0d = g.Sx0'*g.Gd*g.Sx + g.Sy0'*g.Gd*g.Sy;
Q00d = g.Sx0'*g.Gd*g.Sx0 + g.Sy0'*g.Gd*g.Sy0;

% uses acceleration not velocity in effort term
Qd = g.Sx'*g.Gdd*g.Sx + g.Sy'*g.Gdd*g.Sy;
Q0d = g.Sx0'*g.Gdd*g.Sx + g.Sy0'*g.Gdd*g.Sy;
Q00d = g.Sx0'*g.Gdd*g.Sx0 + g.Sy0'*g.Gdd*g.Sy0;

Q1d = g.Sx'*g.G1d*g.Sx + g.Sy'*g.G1d*g.Sy;
Q01d = g.Sx0'*g.G1d*g.Sx + g.Sy0'*g.G1d*g.Sy;

%Ji = 1/2*(ui'*t.Rii*ui + 2*u_i'*t.Ri_i*ui + u_i'*t.R_i_i*u_i) + ...
%   t.ri*ui + t.r_i*u_i + t.z;
%fmax5 %fmax 0.1
Riiterm1_1 = ui'*g.w1*QVP{n}*ui ;%209;206 %727;1345
Riiterm2_1 = ui'*(g.w2 + g.r*g.stiff^2)*Q *ui;%619;589 % 5835;10058
Riiterm3_1= ui'*g.mass^2*g.r*Qd*ui;%0.65;0.78 %2272; 4107
Riiterm4_1= ui'*2*g.mass*g.stiff*g.r*Q1d*ui;%-0.5;-0.77; %-3213;-7660
Riiterm1 = Riiterm1_1 + Riiterm2_1 + Riiterm3_1 + Riiterm4_1

Ri_i_term1_1= -ui'*(g.w2 + g.r*g.stiff^2)*Q*u_i ;%-619;-589 % -5578;-9868
Ri_i_term2_1= -ui'*(g.r*g.mass*g.stiff*Q1d)*u_i;%0.25;0.38 %1382;3210
Ri_i_term1 = Ri_i_term1_1+Ri_i_term2_1

R_i_i_term1_1= u_i'*(g.w2 + g.r*g.stiff^2)*Q*u_i;%%-619;-589 %-5835;-10058
R_i_i_term1 = R_i_i_term1_1

R1 = Riiterm1 + 2*Ri_i_term1 + R_i_i_term1 %210;207 %3067;5638

ri_term1_1= + 0.5*g.w1*(+ g.u0'*QVP0{n})*ui;%-0.8;0.788; -0.81;2.6
ri_term2_1= + 0.5*g.w1*(-2*g.VP1(1)'*g.giVP{n}*g.Sx) *ui; %-64.79;-63 %-65;-209
ri_term3_1= + 0.5*g.w1*(-2*g.VP1(2)'*g.giVP{n}*g.Sy) *ui;%-155.02;-144 %-308;-308
ri_term4_1= g.mass^2*g.r*g.u0'*Q0d*ui;%-0.57;-0.58; %1427;-1534
ri_term5_1= g.r*g.mass*g.stiff*g.u0'*Q01d*ui;

ri_term1 = ri_term1_1+ri_term2_1+ri_term3_1+ri_term4_1+ri_term5_1

r_i_term1_1 = (- g.r*g.mass*g.stiff*g.u0'*Q01d)*u_i; %-0.14;-0.14 %-374;-333
r_i_term1 = r_i_term1_1

r1 = ri_term1 + r_i_term1 % ; % -1801;-2049

z_1 = (0.5*(g.w1*(g.u0'*QVP00{n}*g.u0 ...
    + g.VP1(1).^2 + g.VP1(2).^2 ...
    - 2*g.VP1(1)*g.giVP{n}*g.Sx0*g.u0 - 2*g.VP1(2)*g.giVP{n}*g.Sy0*g.u0) ...
    +g.mass^2*g.r*g.u0'*Q00d*g.u0)) %105;102 %780;777
J1 = 0.5*R1 + r1 + z_1

% Player 2
Riiterm1_2= u_i'*g.w1*QVP{n}*u_i ;
Riiterm2_2= u_i'*(g.w2 + g.r*g.stiff^2)*Q *u_i;
Riiterm3_2= u_i'*g.mass^2*g.r*Qd*u_i;
Riiterm4_2= u_i'*2*g.mass*g.stiff*g.r*Q1d*u_i;
Riiterm2 = Riiterm1_2+Riiterm2_2+Riiterm3_2+Riiterm4_2

Ri_i_term1_2= -u_i'*(g.w2 + g.r*g.stiff^2)*Q*ui ;
Ri_i_term2_2 = -u_i'*(g.r*g.mass*g.stiff*Q1d)*ui;
Ri_i_term2 = Ri_i_term1_2+Ri_i_term2_2

R_i_i_term1_2 = ui'*(g.w2 + g.r*g.stiff^2)*Q*ui;
R_i_i_term2 =  R_i_i_term1_2

R2 = Riiterm2 + 2*Ri_i_term2 + R_i_i_term2


ri_term1_2 = + 0.5*g.w1*(+ g.u0'*QVP0{n})*u_i;
ri_term2_2 = + 0.5*g.w1*(-2*g.VP2(1)'*g.giVP{n}*g.Sx) *u_i;
ri_term3_2= + 0.5*g.w1*(-2*g.VP2(2)'*g.giVP{n}*g.Sy) *u_i;
ri_term4_2 = g.mass^2*g.r*g.u0'*Q0d*u_i;
ri_term5_2= g.r*g.mass*g.stiff*g.u0'*Q01d*u_i;
ri_term2 = ri_term1_2+ri_term2_2+ri_term3_2+ri_term4_2+ri_term5_2

r_i_term1_2= (- g.r*g.mass*g.stiff*g.u0'*Q01d)*ui;
r_i_term2 = r_i_term1_2;

r2 = ri_term2 + r_i_term2 

z_2= (0.5*(g.w1*(g.u0'*QVP00{n}*g.u0 ...
    + g.VP2(1).^2 + g.VP2(2).^2 ...
    - 2*g.VP2(1)*g.giVP{n}*g.Sx0*g.u0 - 2*g.VP2(2)*g.giVP{n}*g.Sy0*g.u0) ...
    +g.mass^2*g.r*g.u0'*Q00d*g.u0))
J2 =  0.5*R2 +r2 + z_2

Jiok = (J1>0) && (J2>0)
