function e = observer(A,B,C,D,SigmaX,SigmaY,x0,P0)
% Estimator Costructor

e.A = A;
e.C = C;
e.B = B;
e.D = D;

e.SigmaX = SigmaX;
e.SigmaY = SigmaY;

e.x0 = x0;
e.P0 = P0;

xsize = length(A);
ysize = size(C,1);
e.xsize = xsize;
e.ysize = ysize;

e = class(e,'observer');

end

