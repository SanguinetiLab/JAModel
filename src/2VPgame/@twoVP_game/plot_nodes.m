function plot_nodes(g,nn,u1,u1hat,u2,u2hat)

%% All nodes
colors = rand(nn*2,3);
figure
set(gcf,'pos',[0 0 500 600])
subplot(411)
hold on
for i = 1:nn*2
    plot(u1(:,i),'color',colors(i,:),'linewidth',1.5)
    plot(u1hat(:,i),':','color',colors(i,:),'linewidth',1.5)
end
title('U1 - Uhat1')

% Visualization Player 2
subplot(412)
hold on
for i = 1:nn*2
    plot(u2(:,i),'color',colors(i,:),'linewidth',1.5)
    plot(u2hat(:,i),':','color',colors(i,:),'linewidth',1.5)
end
title('U2 - Uhat2')


subplot(413)
hold on
for i = 2:2:nn*2
    plot(u1(:,i),'color',colors(i,:),'linewidth',1.5)
    plot(u1hat(:,i),':','color',colors(i,:),'linewidth',1.5)
end
ylabel('Ynodes')
title('U1 - Uhat1')

% Visualization Player 2
subplot(414)
hold on
for i = 2:2:nn*2
    plot(u2(:,i),'color',colors(i,:),'linewidth',1.5)
    plot(u2hat(:,i),':','color',colors(i,:),'linewidth',1.5)
end
title('U2 - Uhat2')
ylabel('Ynodes')
