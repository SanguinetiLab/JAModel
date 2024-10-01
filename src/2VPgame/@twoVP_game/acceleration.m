function [a1x,a1y,a2x,a2y]=acceleration(g,u1,u2)

% Compute acceleration
a1x = g.giddot*(g.Sx*u1+g.Sx0*g.u0);
a1y = g.giddot*(g.Sy*u1+g.Sy0*g.u0);

a2x = g.giddot*(g.Sx*u2+g.Sx0*g.u0);
a2y = g.giddot*(g.Sy*u2+g.Sy0*g.u0);

