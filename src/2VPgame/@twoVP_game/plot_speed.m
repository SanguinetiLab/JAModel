function plot_speed(g,u1,u2)

figure
set(gcf,'pos',[100 100 150 150])
[pd1x,pd1y,pd2x,pd2y]=velocity(g,u1,u2);

ti = (0:g.dt:g.duration)';
for t=1:size(u1,2)
    line(ti,sqrt(pd1x(:,t).^2+pd1y(:,t).^2),'col','b','lines','-')
    line(ti,sqrt(pd2x(:,t).^2+pd2y(:,t).^2),'col','r','lines','-')
end
ylim([0 0.2])
xlim([0 g.duration])
xlabel('time [s]')
ylabel('speed [m/s]')

end

