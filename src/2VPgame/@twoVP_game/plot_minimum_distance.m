function plot_minimum_distance(g,u1,u2)

[p1x,p1y,p2x,p2y]=trajectory(g,u1,u2);

T = size(p1x,2);
for t=1:T
    md11(t) = min(sqrt((p1x(:,t)-g.VP1(1)).^2+ (p1y(:,t)-g.VP1(2)).^2));
    md12(t) = min(sqrt((p1x(:,t)-g.VP2(1)).^2+ (p1y(:,t)-g.VP2(2)).^2));
    
    md21(t) = min(sqrt((p2x(:,t)-g.VP1(1)).^2+ (p2y(:,t)-g.VP1(2)).^2));
    md22(t) = min(sqrt((p2x(:,t)-g.VP2(1)).^2+ (p2y(:,t)-g.VP2(2)).^2));
end

figure
set(gcf,'pos',[300 100 200 200])
subplot(2,2,1)
scatter(1:T,md11,'b.')
ylim([0 0.05])
title('MD_{11}')
set(gca,'fontname','Times')

subplot(2,2,2)
scatter(1:T,md12,'b.')
ylim([0 0.05])
title('MD_{12}')
set(gca,'fontname','Times')

subplot(2,2,3)
scatter(1:T,md21,'r.')
ylim([0 0.05])
title('MD_{21}')
set(gca,'fontname','Times')

subplot(2,2,4)
scatter(1:T,md22,'r.')
ylim([0 0.05])
title('MD_{22}')
set(gca,'fontname','Times')