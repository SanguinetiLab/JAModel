function plot_trajectory(g,u1,u2,u1ne,u2ne)

% Compute trajectory
p1x = g.gi*(g.Sx*u1+g.Sx0*g.u0);
p1y = g.gi*(g.Sy*u1+g.Sy0*g.u0);

p2x = g.gi*(g.Sx*u2+g.Sx0*g.u0);
p2y = g.gi*(g.Sy*u2+g.Sy0*g.u0);

% Plot trajectory
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
     hl(3)=line(u1(1:2:end,t),u1(2:2:end,t),'col','b','lines','none','marker','o');
     hl(4)=line(u2(1:2:end,t),u2(2:2:end,t),'col','r','lines','none','marker','o');
     ht = text(0,0.05,sprintf('trial: %d',t));
    else
     set(hl(1),'xdata',p1x(:,t),'ydata',p1y(:,t));
     set(hl(2),'xdata',p2x(:,t),'ydata',p2y(:,t));
     set(hl(3),'xdata',u1(1:2:end,t),'ydata',u1(2:2:end,t));
     set(hl(4),'xdata',u2(1:2:end,t),'ydata',u2(2:2:end,t));
     
     set(ht,'string',sprintf('trial: %d',t));
    end
    xlim([-0.06 0.06])
    ylim([-0.06 0.06]) 
    pause(0.1)
end

nt = get(get(g,'task1'),'nt');

if nt==1
    
    p1nex = g.gi*(g.Sx*u1ne+g.Sx0*g.u0);
    p1ney = g.gi*(g.Sy*u1ne+g.Sy0*g.u0);
    
    p2nex = g.gi*(g.Sx*u2ne+g.Sx0*g.u0);
    p2ney = g.gi*(g.Sy*u2ne+g.Sy0*g.u0);
    
    line(p1nex,p1ney,'col','b','lines','--');
    line(p2nex,p2ney,'col','r','lines','--');
    line(u1ne(1:2:end),u1ne(2:2:end),'col','b','marker','*','lines','none')
    line(u2ne(1:2:end),u2ne(2:2:end),'col','r','marker','*','lines','none')
else
    
    for  n=1:nt
        p1nex = g.gi*(g.Sx*u1ne{n}+g.Sx0*g.u0);
        p1ney = g.gi*(g.Sy*u1ne{n}+g.Sy0*g.u0);
        
        p2nex = g.gi*(g.Sx*u2ne{n}+g.Sx0*g.u0);
        p2ney = g.gi*(g.Sy*u2ne{n}+g.Sy0*g.u0);
        
        line(p1nex,p1ney,'col','b','lines',':');
        line(p2nex,p2ney,'col','r','lines',':');
        line(u1ne{n}(1:2:end),u1ne{n}(2:2:end),'col','b','marker','*','lines','none')
        line(u2ne{n}(1:2:end),u2ne{n}(2:2:end),'col','r','marker','*','lines','none')
        
    end
end
xlim([-0.06 0.06])
ylim([-0.06 0.06])   
