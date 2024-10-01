function g = twoVP_game(M,VP1,VP2,tcE,tcL,start,final,duration,VPradius,tgtradius)
% DEFINE_2VP_TASK defines the parameters of a 2-VP task

% Define actions
% u1 = [p1(tc1) p1(tc2)]
% u2 = [p2(tc1) p2(tc2)]

tcM = 0.5;

g.M=M; % M is the number of intermediate nodes
g.VP1 = VP1;
g.VP2 = VP2;
g.start = start;
g.final = final;
g.duration = duration;
g.VPradius = VPradius;
g.tgtradius=tgtradius;

g.tcE = tcE*duration;
g.tcL = tcL*duration;
g.tcM = tcM*duration;

g.u0 = [start final]';

% Create spline approximation of trajectory
g.dt = 0.01; % time resolution
g.dt = 0.01; % time resolution


% spline nodes...
% first two and last two to make sure start and end velocities are near zero
nodes_times = linspace(0,g.duration,M+2);
g.times = [0 g.dt nodes_times(2:(M+1)) g.duration-g.dt g.duration]';
%g.times = [0 0.2 nodes_times(2:(M+1)) g.duration-0.2 g.duration]';

g.t = (0:g.dt:g.duration)';
g.K = length(g.t);

phs = phspline(g.times);
g.S = get(phs,'S');

g.gi = gmtx(phs,g.t);
g.G = g.gi'*g.gi;

g.gidot = gdmtx(phs,g.t);
g.Gd = g.gidot'*g.gidot;

g.giddot = gddmtx(phs,g.t);
g.Gdd = g.giddot'*g.giddot;

g.Gdd_1 = g.giddot'*g.gi;

g.Gdd_d = g.giddot'*g.gidot;
g.Gd_dd = g.gidot'*g.giddot; %maybe better for matrix size??
g.G1_d = g.gi'*g.gidot;
g.Gd_1 = g.gidot'*g.gi; %maybe better for matrix size??

% NB1: In the matrices below the rows indicate the nodes, the columns
% indicate:
% in Px, Py the components (x and y) of u1/u2 (the variable nodes)
% in Px0, Py0 the components (x and y) of u0 (the fixed nodes, i.e. start and end)

% In general, there are M+4 nodes where the first two are start and
% the last two are end (each repeated twice).
% M is the number of variable nodes
% In conclusion, all the matrices below have M+4 rows
% Px and Py must have 2M columns, Px0 and Py0 have 4 columns

Px = zeros(M+4,2*M);
Py = zeros(M+4,2*M);
for m = 1:M
    Px(m+2,2*m-1) = 1;
    Py(m+2,2*m) = 1;
end

Px0 = zeros(M+4,4);%2*M);
Px0([1 2],1) = 1; % start_x
Px0([end-1 end],end-1) = 1; % final_x

Py0 = zeros(M+4,4);%2*M);
Py0([1 2],2) = 1; % start_y
Py0([end-1 end],end) = 1; % final_y

g.Sx = g.S*Px;
g.Sy = g.S*Py;
g.Sx0 = g.S*Px0;
g.Sy0 = g.S*Py0;

% This approximation allows to define the task (cost functional) in terms
% of the whole trajectory pix(t), piy(t). We have:
%
% pix(t) = gi(t)*(Sx*ui+Sx0*u0);
% piy(t) = gi(t)*(Sy*ui+Sy0*u0);
%
% ... and its derivatives, ie:
%
% pidx(t) = gid(t)*(Sx*ui+Sx0*u0);
% pidy(t) = gid(t)*(Sy*ui+Sy0*u0);
%
% piddx(t) = gidd(t)*(Sx*ui+Sx0*u0);
% piddy(t) = gidd(t)*(Sy*ui+Sy0*u0);

% NB2: In order to fully generalize to any number of nodes, M, we would
% like to specify the crossing times tc1 and tc2 independently of the
% times at which the nodes are placed (these should be uniformly
% spaced as they determine the frequency content of the trajectory.
% This requires some rearrangements of the constructor and of the 'w1'
% part of the cost function below.

% The position at time tci is defined as:
% pix(tc) = gi(tc)*(Sx*ui+Sx0*u0);
% piy(tc) = gi(tc)*(Sy*ui+Sy0*u0);

% Cost function weights (Bryson's rule)
g.w1 = 1/VPradius.^2;

dmax = 0.00015; % maximum distance between partners
dmax = 0.0015; % 0.015 maximum distance between partners
dmax = 0.002;

connected = 1;
% connected = 0;

g.w2 = connected/(dmax.^2)/g.K;
%g.w2 = 0;
g.visc = .1*connected; %0.1*connected; 
g.stiff = 80*connected; %80*connected;
g.mass = 3; %3
fmax = 1; %1
g.r = 5*(1./(fmax)^2/g.K); %5*

% Player 1:
% J_1(u_1,u_2, s_1) = 1/2 [w_1 (p1x(tcj)-VP(1))^T (p1x(tcj)-VP(1)) +
%                        + w_1 (p1y(tcj)-VP(2))^T (p1y(tcj)-VP(2)) +
%
%                        + w_2 1/K \sum_{k=0}{K} (p1x(k)-px2(k))^T(px1(k)-px2(k)) +
%                        + w_2 1/K \sum_{k=0}{K} (p1y(k)-py2(k))^T(py1(k)-py2(k)) +
%
%                        + r   1/K \sum_{k=0}{K} F1x(k)^T F1x(k) +
%                        + r   1/K \sum_{k=0}{K} F1y(k)^T F1y(k)]

% where F1x(k) = m p1ddx(k) + k (p1x(k)-p2x(k)) + b (p1dx(k)) therefore:

% F1x(k)^T F1x(k) = m^2 * p1ddx(k)^T*p1ddx(k)
%                 + k^2 * (p1x(k)-p2x(k))^T*(p1x(k)-p2x(k))
%                 + b^2 * pidx(k)^T*pidx(k)
%                 + 2mk * p1ddx(k)^T*(p1x(k)-p2x(k))
%                 + 2mb * p1ddx(k)^T*p1dx(k)
%                 + 2kb * (p1x(k)-p2x(k))^T*p1dx(k)   %%%%%ok

% J_1(u_1,u_2, s_1) = 1/2 [w_1 (p1x(tcj)-VP(1))^T (p1x(tcj)-VP(1)) +
%                          w_1 (p1y(tcj)-VP(2))^T (p1y(tcj)-VP(2)) +
%
%                        + w_2/K \sum_{k=0}{K} (p1x(k)-p2x(k))^T(p1x(k)-p2x(k))+
%                        + w_2/K \sum_{k=0}{K} (p1y(k)-p2y(k))^T(p1y(k)-p2y(k))+
%
%                        + r/K \sum_{k=0}{K} ( m^2 p1ddx(k)^T*p1ddx(k) +
%                                              k^2 (p1x(k)-p2x(k))^T*(p1x(k)-p2x(k))+
%                                              b^2 pidx(k)^T*pidx(k) +
%                                              2mk p1ddx(k)^T(p1x(k)-p2x(k)) +
%                                              2mb p1ddx(k)^T*p1dx(k)) +
%                                              2kb (p1x(k)-p2x(k))^T*p1dx(k) ) +
%
%                        + r/K \sum_{k=0}{K} ( m^2 p1ddy(k)^T*p1ddy(k) +
%                                              k^2 (p1y(k)-p2y(k))^T(p1y(k)-p2y(k))+
%                                              b^2 pidy(k)^T*pidy(k) +
%                                              2mk p1ddy(k)^T(p1y(k)-p2y(k))+
%                                              2mb p1ddy(k)^T*p1dy(k) +
%                                              2kb (p1y(k)-p2y(k))^T*p1dy(k) )]
%                                              %%%%%%%%ok

%             = 1/2*w_1*[(Sx*u_1+Sx0*u0)'*G_VPj*(Sx*u_1+Sx0*u0)
%                       +(Sy*u_1+Sy0*u0)'*G_VPj*(Sy*u_1+Sy0*u0)
%                       +VP1(1)^2 + VP1(2)^2
%                       -2*VP1(1) * gVPj*(Sx *u_1 + Sx0* u0)
%                       -2*VP1(2) * gVPj*(Sy *u_1 + Sy0* u0)]
%
%               1/2 * w_2 * (Sx *u_1 - Sx* u_2)' G (Sx u_1 - Sx * u_2) +
%               1/2 * w_2 * (Sy *u_1 - Sy* u_2)' G (Sy u_1 - Sy * u_2) +
%
%               1/2 * r* m^2 *(Sx *u_1 + Sx0* u0)' Gdd (Sx *u_1 + Sx0* u0) +
%               1/2 * r* m^2 *(Sy *u_1 + Sy0* u0)' Gdd (Sy *u_1 + Sy0* u0) +
%
%               1/2 * r* k^2 *(Sx *u_1 - Sx* u_2)' G (Sx u_1 - Sx * u_2) +
%               1/2 * r* k^2 *(Sy *u_1 - Sy* u_2)' G (Sy u_1 - Sy * u_2) +
%
%               1/2 * r* b^2 *(Sx *u_1 + Sx0* u0)' Gd (Sx u_1 + Sx0* u0) +
%               1/2 * r* b^2 *(Sy *u_1 + Sy0* u0)' Gd (Sy u_1 + Sy0* u0) +
%
%               1/2 * r* 2*m*k * (Sx *u_1 + Sx0* u0)' Gdd_1 (Sx u_1 - Sx * u_2) +
%               1/2 * r* 2*m*k * (Sy *u_1 + Sy0* u0)' Gdd_1 (Sy u_1 - Sy * u_2) +
%
%               1/2 * r* 2*m*b * (Sx *u_1 + Sx0* u0)' Gdd_d (Sx *u_1 + Sx0* u0) +
%               1/2 * r* 2*m*b * (Sy *u_1 + Sy0* u0)' Gdd_d (Sy *u_1 + Sy0* u0) +
%
%               1/2 * r* 2*k*b * (Sx *u_1 - Sx* u_2)' G1_d (Sx *u_1 + Sx0* u0) +
%               1/2 * r* 2*k*b * (Sy *u_1 - Sx* u_2)' G1_d (Sy *u_1 + Sy0* u0) = %ok

%             = 1/2* u1' * [w_1         * (Sx' G_VPj Sx + Sy' G_VPj Sy)
%                          (w_2+rk^2)   * (Sx' G Sx + Sy' G Sy)+
%                           rm^2        * (Sx'Gdd Sx + Sy' Gdd Sy)+
%                           rb^2        * (Sx'Gd Sx + Sy' Gd Sy) + 
%                           2mkr        * (Sx' Gdd_1 Sx + Sy' Gdd_1 Sy) + 
%                           2mbr        * (Sx' Gdd_d Sx + Sy' Gdd_d Sy) + 
%                           2kbr        * (Sx' G1_d Sx + Sy' G1_d Sy)] * u1 +
%
%               1/2* u1' *[ -2(w_2+rk^2)* (Sx'G Sx + Sy'G Sy) +
%                           -2rmk       * (Sx'Gdd_1 Sx + Sy'Gdd_1 Sy) + 
%                           -2rkb       * (Sx'G1_d Sx + Sy'G1_d Sy)] * u2+
%
%               1/2* u2' *[(w_2+rk^2)   * (Sx' G Sx + Sy' G Sy)] * u2 +
%
%                    { w_1         * [ u0'*(Sx0' GVPj Sx + Sy0' GVPj Sy) +
%                                      -VP1(1)'*gVPj*Sx - VP1(2)'*gVPj*Sy] +
%                     rm^2         * u0'*(Sx0' Gdd Sx + Sy0' Gdd Sy) +
%                     rb^2         * u0'*(Sx0' Gd Sx + Sy0' Gd Sy) +
%                     rmk          * u0'*(Sx0' Gdd_1 Sx + Sy0' Gdd_1 Sy) + 
%                     2rmb          * u0'*(Sx0' Gdd_d Sx + Sy0' Gdd_d Sy) + 
%                     rkb          * u0'*(Sx0' G1_d Sx + Sy0' G1_d Sy)
%                     } * u1 +
%
%                    {-rmk         * u0'*(Sx0' Gdd_1 Sx + Sy0' Gdd_1 Sy) + 
%                     -rbk         * u0'*(Sx0' G1_d Sx + Sy0' G1_d Sy)} * u2 +
%
%               1/2 * w_1           * [ VP1(1)^2+VP1(2)^2 +
%                                       u0'*(Sx0' GVPj Sx0 + Sy0' GVPj Sy0)*u0 +
%                                     - 2*[VP1(1) gVPj Sx0 + VP1(2) gVPj Sy0]*u0 ] +
%               1/2 * rm^2          * u0' *(Sx0' Gdd Sx0 + Sy0' Gdd Sy0) u0 + 
%               1/2 * rb^2          * u0' *(Sx0' Gd Sx0 + Sy0' Gd Sy0) u0 +
%               1/2 * 2rmb          * u0' *(Sx0' Gdd_d Sx0 + Sy0' Gdd_d Sy0) * u0 = %ok

% where Qj may be QE (crossing early) or QL (crossing late):


% for VP part
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
Qdd = g.Sx'*g.Gdd*g.Sx + g.Sy'*g.Gdd*g.Sy;
Q0dd = g.Sx0'*g.Gdd*g.Sx + g.Sy0'*g.Gdd*g.Sy;
Q00dd = g.Sx0'*g.Gdd*g.Sx0 + g.Sy0'*g.Gdd*g.Sy0;

Qdd_1 = g.Sx'*g.Gdd_1*g.Sx + g.Sy'*g.Gdd_1*g.Sy;
Q0dd_1 = g.Sx0'*g.Gdd_1*g.Sx + g.Sy0'*g.Gdd_1*g.Sy;

Qdd_d = g.Sx'*g.Gdd_d*g.Sx + g.Sy'*g.Gdd_d*g.Sy;
Q0dd_d = g.Sx0'*g.Gdd_d*g.Sx + g.Sy0'*g.Gdd_d*g.Sy;
Q00dd_d = g.Sx0'*g.Gdd_d*g.Sx0 + g.Sy0'*g.Gdd_d*g.Sy0;

Q1_d = g.Sx'*g.G1_d*g.Sx + g.Sy'*g.G1_d*g.Sy;
Q01_d = g.Sx0'*g.G1_d*g.Sx + g.Sy0'*g.G1_d*g.Sy;
Q001_d = g.Sx0'*g.G1_d*g.Sx0 + g.Sy0'*g.G1_d*g.Sy0;

% J_1(u_1,u_2, s_1) =
%             = 1/2 *[u1'* [w_1*QVPj + (w_2+rk^2)*Q + rm^2*Qdd + rb^2*Qd + 2mkr*Qdd_1 + 2mbr*Qdd_d + 2kbr*Q1_d] *u1 + ok
%
%                      - 2 u1'* [(w_2+rk^2)*Q + rmk*Qdd_1 + rkb*Q1_d] *u2 
%
%                      + u2'* (w_2+rk^2)*Q *u2 ] + ok
%
%               +[w_1 * ( 0.5*u0' Q0VPj +
%                               - VP1(1) g.giE Sx - VP1(2) g.giE Sy) +
%                        rm^2 * u0'*Q0dd +
%                        rb^2 * u0'*Q0d + 
%                        rmk * u0'*Q0dd_1 + 
%                        2rmb * u0'*Q0dd_d + 
%                        rkb * u0'*Q01_d] u1 +
%
%                - (rmk * u0'*Qdd_1 + rbk*u0'*Q01_d)*u2 +
%
%               1/2 * {w1   * [ VP(1)^2 + VP(2)^2 +
%                               u0' QVP00 u0 +
%                               - 2 (VP(1)*g.giE*Sx0 - 2VP(2)*g.giE*Sy0 )u0  ]+
%                    rm^2    *  u0' Q00dd u0 + 
%                    rb^2    *  u0' Q00d u0 + 
%                    2rmb    *  u0' Q00dd_d u0} % ok

% where Qj may be QE (crossing early) or QL (crossing late):

s2 = [2 1 3]; % order of crossing times in player 2: late early center
%s2 = [1 2 3];

for s = 1:3
    % Player 1: strategy 1: VP1 crossed early;
    %           strategy 2: VP1 crossed late
    %           strategy 3: VP1 crossed middle
    Rii_1{s} =  g.w1*QVP{s}...
        + (g.w2  + g.r*g.stiff^2)*Q...
        +  g.r*g.mass^2*Qdd...
        +  g.r*g.visc^2*Qd...
        +  2*g.r*g.mass*g.stiff*Qdd_1...%+
        +  2*g.r*g.mass*g.visc*Qdd_d... %+
        +  2*g.r*g.stiff*g.visc*Q1_d; %+
 
    Ri_i_1{s} = - (g.w2+g.r*g.stiff^2)*Q...
        -  g.r*g.mass*g.stiff*Qdd_1...%-
        -  g.r*g.stiff*g.visc*Q1_d; %-

    R_i_i_1{s} = (g.w2 + g.r*g.stiff^2)*Q;

    ri_1{s} = g.w1*(g.u0'*QVP0{s}...
        - g.VP1(1)'*g.giVP{s}*g.Sx - g.VP1(2)'*g.giVP{s}*g.Sy)...
        +g.r*g.mass^2*g.u0'*Q0dd...
        +g.r*g.visc^2*g.u0'*Q0d...
        +g.r*g.mass*g.stiff*g.u0'*Q0dd_1...%+
        +2*g.r*g.mass*g.visc*g.u0'*Q0dd_d...%+
        +g.r*g.stiff*g.visc*g.u0'*Q01_d;%+


    r_i_1{s} = -g.r*g.mass*g.stiff*g.u0'*Q0dd_1...%-
        -g.r*g.stiff*g.visc*g.u0'*Q01_d;%-

    z_1{s} = 0.5*g.w1*(g.u0'*QVP00{s}*g.u0 ...
                        + g.VP1(1).^2 + g.VP1(2).^2 ...
                        - 2*VP1(1)*g.giVP{s}*g.Sx0*g.u0 ...
                        - 2*VP1(2)*g.giVP{s}*g.Sy0*g.u0) ...
            +0.5*g.r*g.mass^2*g.u0'*Q00dd*g.u0...
            +0.5*g.r*g.visc^2*g.u0'*Q00d*g.u0...
            +g.r*g.mass*g.visc*g.u0'*Q00dd_d*g.u0;%+

    % Player 2: strategy 1: VP2 crossed late;
    %           strategy 2: VP2 crossed early;
    %           strategy 3: VP2 crossed middle;

    Rii_2{(s)} =  g.w1*QVP{s2(s)}...
        + (g.w2 + g.r*g.stiff^2)*Q...
        +  g.r*g.mass^2*Qdd...
        +  g.r*g.visc^2*Qd...
        +  2*g.r*g.mass*g.stiff*Qdd_1...%+
        +  2*g.r*g.mass*g.visc*Qdd_d...%+
        +  2*g.r*g.stiff*g.visc*Q1_d;%+

    Ri_i_2{(s)} = - (g.w2+g.r*g.stiff^2)*Q...
        -  g.r*g.mass*g.stiff*Qdd_1...%-
        -  g.r*g.stiff*g.visc*Q1_d;%-

    R_i_i_2{(s)} = (g.w2 + g.r*g.stiff^2)*Q;

    ri_2{(s)} = g.w1*(g.u0'*QVP0{s2(s)}...
        -g.VP2(1)'*g.giVP{s2(s)}*g.Sx...
        -g.VP2(2)'*g.giVP{s2(s)}*g.Sy)...
        +g.r*g.mass^2*g.u0'*Q0dd...
        +g.r*g.visc^2*g.u0'*Q0d...
        +g.r*g.mass*g.stiff*g.u0'*Q0dd_1...%+
        +2*g.r*g.mass*g.visc*g.u0'*Q0dd_d...%+
        +g.r*g.stiff*g.visc*g.u0'*Q01_d;%+

    r_i_2{(s)} = -g.r*g.mass*g.stiff*g.u0'*Q0dd_1...%-
        -g.r*g.stiff*g.visc*g.u0'*Q01_d;%-

    z_2{(s)} = 0.5*g.w1*(g.u0'*QVP00{s2(s)}*g.u0 ...
        + g.VP2(1).^2 + g.VP2(2).^2 ...
        - 2*VP2(1)*g.giVP{s2(s)}*g.Sx0*g.u0 ...
        - 2*VP2(2)*g.giVP{s2(s)}*g.Sy0*g.u0) ...
        +0.5*g.r*g.mass^2*g.u0'*Q00dd*g.u0...
        +0.5*g.r*g.visc^2*g.u0'*Q00d*g.u0...
        +g.r*g.mass*g.visc*g.u0'*Q00dd_d*g.u0;%+

end

g.task1 = task(Rii_1,Ri_i_1,R_i_i_1,ri_1,r_i_1,z_1);
g.task2 = task(Rii_2,Ri_i_2,R_i_i_2,ri_2,r_i_2,z_2);

% % No decision - GNB
% g.task1 = task(Rii_1{1},Ri_i_1{1},R_i_i_1{1},ri_1{1},r_i_1{1},z_1{1});
% g.task2 = task(Rii_2{1},Ri_i_2{1},R_i_i_2{1},ri_2{1},r_i_2{1},z_2{1});

g = class(g,'twoVP_game');


