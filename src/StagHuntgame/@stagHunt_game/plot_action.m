function plot_action(g,traj1,traj2,sigcol,action_range)


norm_factor = max([sigcol]);
sigcol_new = sigcol./norm_factor ;

figure
set(gcf,'pos',[0 0 500 250])

% subplot(131)
[sorted_sigcol, sortidx_sigcol] = sort(sigcol_new);
xs1 = traj1;%(sortidx_sigcol);
ys1 = traj2;%(sortidx_sigcol);
cmap = jet(length(xs1));

scatter(xs1, ys1, 20, [0,0,128]./255, 'filled')
scatter(xs1, ys1, 20, cmap, 'filled')

xlim([action_range(1) action_range(end)])
%ylim([-action_range(end) -action_range(1) ])
ylim([action_range(1) action_range(end)])

xticks([-1 0 1])
yticks([-1 0 1])
%yticks([-16 -12 -8 -4 -1])
%yticklabels({'16','12','8','4','1'})
xlabel('$u_1$','fontname','times','interpreter','latex')
ylabel('$u_2$','fontname','times','interpreter','latex')

set(gca,'fontname','times')
set(gcf,'pos',[0 0 150 150])
axis square
