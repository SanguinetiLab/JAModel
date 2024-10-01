function [d1,d2] = get_traj_differences(g,nn,u1,u1hat,u2,u2hat)


T = size(u1,2);

for t=1:T
    for n = 1:nn
        dist1(nn,t) = mean(sqrt((u1((nn*2-1),t)-u1hat((nn*2-1),t)).^2+ (u1((nn*2),t)-u1hat((nn*2),t)).^2));
        dist2(nn,t) = mean(sqrt((u2((nn*2-1),t)-u2hat((nn*2-1),t)).^2+ (u2((nn*2),t)-u2hat((nn*2),t)).^2));
    end
end

d1 = mean(dist1);
d2 = mean(dist2);