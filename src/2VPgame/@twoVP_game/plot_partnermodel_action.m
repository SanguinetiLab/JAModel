function plot_partnermodel_action(g,x1,x2,u1,u2)
% x1 is partner model of player 1 i.e. player's 2 estimted action
% x2 is partner model of player 1 i.e. player's 1 estimted action
% u1 is partner 1 action
% u2 is partner 2 action

% Partner Model
[px1x,px1y,px2x,px2y]=trajectory(g,x1,x2);

% Action
[pu1x,pu1y,pu2x,pu2y]=trajectory(g,u1,u2);

figure
set(gcf,'pos',[500 100 300 300])

% Workspace
ths = 0:pi/10:(2*pi);
xx = [cos(ths)' sin(ths)'];

xlabel('x (m)')
ylabel('y (m)')

patch(g.VP1(1)+g.VPradius*xx(:,1),g.VP1(2)+g.VPradius*xx(:,2),[0.5 0.5 1])
hold on
patch(g.VP2(1)+g.VPradius*xx(:,1),g.VP2(2)+g.VPradius*xx(:,2),[1 0.5 0.5])

patch(g.start(1)+g.tgtradius*xx(:,1),g.start(2)+g.tgtradius*xx(:,2),'y')
patch(g.final(1)+g.tgtradius*xx(:,1),g.final(2)+g.tgtradius*xx(:,2),'y')

nodes_sz = size(u1,1);


for t = 1:size(u1,2)
    if t ==1
        % Partner Model trajectory
        hl(1)=line(px1x(:,t),px1y(:,t),'col','m','lines',':');
        hl(2)=line(px2x(:,t),px2y(:,t),'col','c','lines',':');
        hl(3)=line(x1([1:2:(nodes_sz-1)],t),x1([2:2:nodes_sz],t),'col','m','lines','none','marker','s');
        hl(4)=line(x2([1:2:(nodes_sz-1)],t),x2([2:2:nodes_sz],t),'col','c','lines','none','marker','s');

        % Actions
        hl(5)=line(pu1x(:,t),pu1y(:,t),'col','b','lines','-');
        hl(6)=line(pu2x(:,t),pu2y(:,t),'col','r','lines','-');
        hl(7)=line(u1([1:2:(nodes_sz-1)]),u1([2:2:nodes_sz]),'col','b','lines','none','marker','o');
        hl(8)=line(u2([1:2:(nodes_sz-1)]),u2([2:2:nodes_sz]),'col','r','lines','none','marker','o');

        ht = text(0,0.05,sprintf('trial: %d',t));
    else
        set(hl(1),'xdata',px1x(:,t),'ydata',px1y(:,t));
        set(hl(2),'xdata',px2x(:,t),'ydata',px2y(:,t));
        set(hl(3),'xdata',x1([1:2:(nodes_sz-1)],t),'ydata',x1([2:2:nodes_sz],t));
        set(hl(4),'xdata',x2([1:2:(nodes_sz-1)],t),'ydata',x2([2:2:nodes_sz],t));

        set(hl(5),'xdata',pu1x(:,t),'ydata',pu1y(:,t));
        set(hl(6),'xdata',pu2x(:,t),'ydata',pu2y(:,t));
        set(hl(7),'xdata',u1([1:2:(nodes_sz-1)],t),'ydata',u1([2:2:nodes_sz],t));
        set(hl(8),'xdata',u2([1:2:(nodes_sz-1)],t),'ydata',u2([2:2:nodes_sz],t));

        set(ht,'string',sprintf('trial: %d',t));
    end
    xlim([-0.06 0.06])
    ylim([-0.06 0.06]) 
    pause(0.1)
end
xlim([-0.06 0.06])
ylim([-0.06 0.06])
