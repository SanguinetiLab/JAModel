function giddot = gddmtx(phs,t)
% GDDMTX  creates the GDDMATRIX gij = [6*(ti-tj) 0 0]  [Tx(M+2)]
% where T = length(t)

T = length(t);
giddot = [6*(t*ones(1,phs.M)-phs.t(:,ones(T,1))').*sign(t*ones(1,phs.M)-phs.t(:,ones(T,1))') zeros(T,1) zeros(T,1)];
    

