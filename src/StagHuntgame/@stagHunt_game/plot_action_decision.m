function plot_action_decision(g,u1,u2,decision1,decision2,action_range)



figure
subplot(211)
plot(u1,'b','linewidth',1.5)
ylabel('Actions','fontname','times','interpreter','latex')
%title(['$u_1$'],'interpreter','latex')
%title('Actions Player1','fontname','times')
%xlabel('Trials','fontname','times')
ylim([action_range(1)-0.1 action_range(end)+0.1])
yline(action_range(end))
yline(action_range(1))
yticks([-1 1])
% yticks([-1 1])
yticklabels({'R','S'})
% box off

%set(gca,'fontname','times')

subplot(212)
plot(u2,'r','linewidth',1.5)
ylabel('Actions')%,'fontname','times','interpreter','latex')
%title(['$u_2$'],'interpreter','latex')
%title('Actions Player2','fontname','times')
xlabel('Trials');%,'fontname','times')
ylim([action_range(1)-0.1 action_range(end)+0.1])
yline(action_range(end))
yline(action_range(1))
yticks([-1 1])
yticklabels({'R','S'})
% box off
%set(gca,'fontname','times')
set(gcf,'pos',[300 0 350 180])

figure
subplot(211)
plot(decision1,'b','linewidth',1.5)
yticks([1 2])
ylim([0.9 2.1])
%title('Decisions Player1','fontname','times')
%xlabel('Trials','fontname','times')
yticklabels({'R','S'})
% box off
% set(gca,'fontname','times')

subplot(212)
plot(decision2,'r','linewidth',1.5)
yticks([1 2])
ylim([0.9 2.1])
title('Decisions Player2','fontname','times')
xlabel('Trials','fontname','times')
yticklabels({'R','S'})
% box off
% set(gca,'fontname','times')
set(gcf,'pos',[300 0 350 180])