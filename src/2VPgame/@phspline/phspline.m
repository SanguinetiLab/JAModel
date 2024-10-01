function phs = phspline(times)
% PHSPLINE defines a poly-harmonic cubic spline 
% PHSPLINE creates a poly-harmonic cubic spline with knots in TIMES:
%
% s(t) = \sum_{i=1}^M w_i |t-times(i)|^3 + w_{M+1} t + W_{M+2}

% 

phs.t = times;
phs.M = length(times);

tt = phs.t(:,ones(1,phs.M));
A = abs(tt-tt').^3;
% for i=1:phs.M
%     for j=1:phs.M
%         A(i,j) = abs(phs.t(i)*ones(1,phs.M)-phs.t(j))^3;
%     end
% end
B = [phs.t ones(phs.M,1)];

T = [A B; B' zeros(2)];

phs.S = inv(T'*T)*T';
phs.S = phs.S(:,1:phs.M);

phs = class(phs,'phspline');

end
