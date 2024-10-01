function gi = gmtx(phs,t)
% GMTX  creates the GMATRIX gij = [abs(t_i-times_j)^3 t_i 1]  [Tx(M+2)]
% where T = length(t)

T = length(t);
gi = [abs(t*ones(1,phs.M)-phs.t(:,ones(T,1))').^3  t ones(T,1)];
    

