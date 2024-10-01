function [K_stat] = StationaryK(e)


% P_stat = A*P_stat*A' + SigmaX - A*P_stat*H'*inv(H*P_stat*H' + SigmaY)*H*P_stat*A'; % Riccati Eq to be solved
% K_stat = A*P_stat*H'*inv(H*P_stat*H' + SigmaY)
A = e.A;
R = e.SigmaY;
S = [];
Q = e.SigmaX;
E = eye(size(A));
C = e.C;

[P_stat, K_stat, eigg] = idare(A,C,Q,R,S,E);
end

