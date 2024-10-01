function plot_animation(g,u1,u2)

% p1x = g.gi*g.S*(g.Px*u1+g.Px0*g.u0);
% p1y = g.gi*g.S*(g.Py*u1+g.Py0*g.u0);
% 
% p2x = g.gi*g.S*(g.Px*u2+g.Px0*g.u0);
% p2y = g.gi*g.S*(g.Py*u2+g.Py0*g.u0);

[p1x,p1y,p2x,p2y]=trajectory(g,u1,u2);


% trajectory
figure
set(gcf,'pos',[100 100 300 300])

ths = 0:pi/10:(2*pi);
xx = [cos(ths)' sin(ths)'];

xlabel('x (m)')
ylabel('y (m)')

patch(g.VP1(1)+g.VPradius*xx(:,1),g.VP1(2)+g.VPradius*xx(:,2),[0.5 0.5 1])
hold on
patch(g.VP2(1)+g.VPradius*xx(:,1),g.VP2(2)+g.VPradius*xx(:,2),[1 0.5 0.5])

patch(g.start(1)+g.tgtradius*xx(:,1),g.start(2)+g.tgtradius*xx(:,2),'y')
patch(g.final(1)+g.tgtradius*xx(:,1),g.final(2)+g.tgtradius*xx(:,2),'y')

for t=1:length(u1)
    if t==1
     hl(1)=line(p1x(:,t),p1y(:,t),'col','b','lines','-');
     hl(2)=line(p2x(:,t),p2y(:,t),'col','r','lines','-');
     hl(3)=line(u1([1 3],t),u1([2 4],t),'col','b','lines','none','marker','o');
     hl(4)=line(u2([1 3],t),u2([2 4],t),'col','r','lines','none','marker','o');
     ht = text(0,0.05,sprintf('trial: %d',t));
    else
     set(hl(1),'xdata',p1x(:,t),'ydata',p1y(:,t));
     set(hl(2),'xdata',p2x(:,t),'ydata',p2y(:,t));
     set(hl(3),'xdata',u1([1 3],t),'ydata',u1([2 4],t));
     set(hl(4),'xdata',u2([1 3],t),'ydata',u2([2 4],t));
     
     set(ht,'string',sprintf('trial: %d',t));
    end
    xlim([-0.01 0.11])
    ylim([-0.06 0.06])
    pause(0.1)
end


xlim([-0.01 0.11])
ylim([-0.06 0.06])   
