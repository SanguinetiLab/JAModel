function plot_md_exp_sim(g,traj1,traj2,sigcol,action_range,tits)

% if sum(diff(sigcol)>0)
norm_factor = max([sigcol]);
sigcol_new = sigcol./norm_factor ;
[sorted_sigcol, sortidx_sigcol] = sort(sigcol_new);

figure
set(gcf,'pos',[700 100 200 200])

for m = 1:size(traj1,1)
    subplot(2,2,m)
    hold on
    plot(action_range,action_range,'-k')
    xs1 = traj1(m,sortidx_sigcol);
    ys1 = traj2(m,sortidx_sigcol);
    cmap = jet(length(xs1));
    scatter(xs1, ys1, 20, cmap, 'filled')
    xlabel('Experiment')
    ylabel('Model')
   % title(tits{m})
    xlim([action_range(1) action_range(end)])
    ylim([action_range(1) action_range(end)])
    axis square
    set(gca,'fontname','Times')
end

figure
set(gcf,'pos',[500 100 200 200])

for m = 1:size(traj1,1)
    if m<=2
        col = 'b'
    else
        col = 'r'
    end
    subplot(2,2,m)
    hold on
    plot(traj1(m,:),'.','Color',col)
%     xlabel('Trials')
%     ylabel('distance')
    title([tits{m} ' Exp'])
    %xlim([action_range(1) action_range(end)])
    ylim([action_range(1) action_range(end)])
    axis square
    set(gca,'fontname','Times')
end

