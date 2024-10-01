function [p1x,p1y,p2x,p2y]=trajectory(g,u1,u2)

% Compute trajectory
p1x = g.gi*(g.Sx*u1+g.Sx0*g.u0);
p1y = g.gi*(g.Sy*u1+g.Sy0*g.u0);

p2x = g.gi*(g.Sx*u2+g.Sx0*g.u0);
p2y = g.gi*(g.Sy*u2+g.Sy0*g.u0);

