function [v1x,v1y,v2x,v2y]=velocity(g,u1,u2)

% Compute velocity
v1x = g.gidot*(g.Sx*u1+g.Sx0*g.u0);
v1y = g.gidot*(g.Sy*u1+g.Sy0*g.u0);

v2x = g.gidot*(g.Sx*u2+g.Sx0*g.u0);
v2y = g.gidot*(g.Sy*u2+g.Sy0*g.u0);

