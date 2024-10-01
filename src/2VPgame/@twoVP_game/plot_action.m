function plot_action(g,u1,u2)
%figure
% Action
[pu1x,pu1y,pu2x,pu2y]=trajectory(g,u1,u2);

ths = 0:pi/10:(2*pi);
xx = [cos(ths)' sin(ths)'];

patch(g.VP1(1)+g.VPradius*xx(:,1),g.VP1(2)+g.VPradius*xx(:,2),'b')
hold on
patch(g.VP2(1)+g.VPradius*xx(:,1),g.VP2(2)+g.VPradius*xx(:,2),'r')
patch(g.start(1)+g.tgtradius*xx(:,1),g.start(2)+g.tgtradius*xx(:,2),'y')
patch(g.final(1)+g.tgtradius*xx(:,1),g.final(2)+g.tgtradius*xx(:,2),'y')

for t=1:size(u1,2)
    % line(pu1x(:,t),pu1y(:,t),'col','b','lines','none','marker','.','markersize',3);
    % line(pu2x(:,t),pu2y(:,t),'col','r','lines','none','marker','.','markersize',3);
    line(pu1x(:,t),pu1y(:,t),'col','b','linewidth',1.5);
    line(pu2x(:,t),pu2y(:,t),'col','r','linewidth',1.5);
    
    line(u1(1:2:end,t),u1(2:2:end,t),'col','b','lines','none','marker','.')
    line(u2(1:2:end,t),u2(2:2:end,t),'col','r','lines','none','marker','.')
end
xlim([-0.07 0.07])%([-0.02 0.12])
ylim([-0.07 0.07])
xlabel('x (m)')
ylabel('y (m)')

axis off


