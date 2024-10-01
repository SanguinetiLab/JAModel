function plot_acceleration(g,u1,u2)

[pdd1x,pdd1y,pdd2x,pdd2y]=acceleration(g,u1,u2);

figure
set(gcf,'pos',[100 100 150 150])
ti = (0:g.dt:g.duration)';
for t=1:size(u1,2)
    line(ti,sqrt(pdd1x(:,t).^2+pdd1y(:,t).^2),'col','b','lines','-')
    line(ti,sqrt(pdd2x(:,t).^2+pdd2y(:,t).^2),'col','r','lines','-')
end
%ylim([0 0.2])
xlim([0 g.duration])
xlabel('time [s]')
ylabel('acceleration [m/s^2]')

end

