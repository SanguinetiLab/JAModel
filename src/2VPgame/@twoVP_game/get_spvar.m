function [spvar1,spvar2] = get_spvar(g,u1,u2,b_trials)
tot_trials = length(u1);
n_epochs = tot_trials/b_trials;

% Spatial Variability from Todorov 2007
for ep = 1:n_epochs
    ep_trials = (1:b_trials) + (ep-1)*b_trials;
    % Undersampling of each trajectory of the targetset
    for t = ep_trials
        [p1x(t,:),p1y(t,:),p2x(t,:),p2y(t,:)]=trajectory(g,u1(:,t),u2(:,t));
    end

    % Mean trajectory computation
    av_tr1_x = mean(p1x,1);
    av_tr1_y = mean(p1y,1);
    av_tr2_x = mean(p2x,1);
    av_tr2_y = mean(p2y,1);

    for i = 1:length(av_tr1_x)
        C1 = cov([p1x(:,i) p1y(:,i)]);
        C2 = cov([p2x(:,i) p2y(:,i)]);
        ell1 = ellipdisp(C1,2);
        var1(i) = (sqrt(trace(C1)))/2;
        ell2 = ellipdisp(C2,2);
        var2(i) = (sqrt(trace(C2)))/2;
    end

    spvar1(:,ep) = mean(var1);
    spvar2(:,ep) = mean(var2);
end

