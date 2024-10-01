function [md11,md12,md21,md22,tc11,tc12,tc21,tc22] = get_minimum_distance(g,u1,u2)

[p1x,p1y,p2x,p2y]=trajectory(g,u1,u2);

T = size(p1x,2);

for t=1:T
    [md11(t),tc11(t)] = min(sqrt((p1x(:,t)-g.VP1(1)).^2+ (p1y(:,t)-g.VP1(2)).^2));
    [md12(t),tc12(t)] = min(sqrt((p1x(:,t)-g.VP2(1)).^2+ (p1y(:,t)-g.VP2(2)).^2));
    
    [md21(t),tc21(t)] = min(sqrt((p2x(:,t)-g.VP1(1)).^2+ (p2y(:,t)-g.VP1(2)).^2));
    [md22(t),tc22(t)] = min(sqrt((p2x(:,t)-g.VP2(1)).^2+ (p2y(:,t)-g.VP2(2)).^2));
end