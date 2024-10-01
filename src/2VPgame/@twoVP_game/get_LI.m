function [LI_11 LI_12 LI_21 LI_22] = get_LI(g,u1,u2,t01,t02,duration)

[p1x,p1y,p2x,p2y]=trajectory(g,u1,u2);
[v1x,v1y,v2x,v2y] = velocity(g,u1,u2);

T = size(p1x,2);

for t=1:T
    tr1 = [p1x(:,t) p1y(:,t)];
    tr2 = [p2x(:,t) p2y(:,t)];
    vel1 = [v1x(:,t) v1y(:,t)];
    vel2 = [v2x(:,t) v2y(:,t)];

    k_spring =150;
    k_visc = 0;

    force1 = -k_spring*(tr1-tr2)-k_visc*vel1;
    force2 = -k_spring*(tr2-tr1)-k_visc*vel2;

    power1 = force1(:,1).*vel1(:,1) + force1(:,2).*vel1(:,2);
    power2 = force2(:,1).*vel2(:,1) + force2(:,2).*vel2(:,2);

    time = linspace(0,duration,size(p1x,1));

    int1 = find(time>=t01-0.300 & time <= t01);
    LI_11(:,t) = mean(power1(int1));
    LI_21(:,t) = mean(power2(int1));

    int2 = find(time>=t02-0.300 & time <= t02);
    LI_12(:,t) = mean(power1(int2));
    LI_22(:,t) = mean(power2(int2));

end