function gidot = gdmtx(phs,t)
% GDMTX  creates the GDMATRIX gij = [3*(ti-tj)^2*sign(ti-tj) 1 0]  [Tx(M+2)]
% where T = length(t)

T = length(t);
gidot = [3*(t*ones(1,phs.M)-phs.t(:,ones(T,1))').^2.*sign(t*ones(1,phs.M)-phs.t(:,ones(T,1))') ones(T,1) zeros(T,1)];
    

