%% Stag Hunt from Yoshida 2009:
clear classes
close all
clear all
clc

addpath(fullfile(cd,'utils'))
figdir = fullfile(cd,'figs');
resdir = fullfile(cd,'results');
check_single_steps = 0;

%% Color specs
purple = [138,43,226]./255;
orange = [230, 149, 0]./255;
light_green = [0.1 0.9 0.1];
green = [0.1 0.6 0.1];
alpha_shade = 0.3;
lw = 1.5;

n_sim = length(dir(resdir))-2;
n_sub = 1;

sigmay_other = 1;%2;%0.1;
n_groups = length(sigmay_other);
sigmax = .2;%2;%0.1;
A = 0.99.*ones(1,n_groups);
lambda = .2.*ones(1,n_groups);%2.*ones(1,n_groups);
lambda_decay = 0.999.*ones(1,n_groups);
P0 = (2*eps.*ones(1,n_groups));
X0 = 0*ones(1,n_groups);


%% Define protocol
n_blocks = 6;
n_reps = 5;
T  = n_blocks*n_reps;
angles = linspace(0, 10*pi,T);
u1imposed = mvnrnd(0,1,T)';%cos(angles);%zeros(n_blocks,n_reps);
u2imposed = u1imposed;%mvnrnd(0,1,T)';%cos(angles);%zeros(n_blocks,n_reps);

for g = 1:n_groups
    for d = 1:n_sub
        %% DEFINE JA MODEL
        %% Define cost functions parameters
        action_range = -2:0.05:2%[0:0.01:1];
        uR = -1;%-1; % Rabbit location
        uS = 1; % Stag location

        u01 = 0;%mean(action_range); % Initial position 1
        u02 = 0;%mean(action_range); %Initial position 2

        % % Gain term
        w = 9/2;
        wR = w; wS =w;
        zS = 1;
        zR = 5;

        sh_game = stagHunt_game(uR, uS, zR, zS, u01,u02,wR,wS);
        task1 = get(sh_game,'task1')
        task2 = get(sh_game,'task2')

        xsz1 = get(task1,'xsize');
        xsz2 = get(task2,'xsize');

        % DISPLAY Nash equilibria
        game = quadraticgame(task1,task2);
        [u1ne,u2ne] = nash_equilibrium(game);

        for n = 1:length(u1ne)
            n
            %plot_action(sh_game,u1ne{n,n},u2ne{n,n},1,action_range)
            J1 = cost(task1,u1ne{n,n},u2ne{n,n});
            J2 = cost(task2,u2ne{n,n},u1ne{n,n});
            %     title(sprintf('NE %d; J1 = %f; J2 = %f',n,J1{n},J2{n}));
        end
        clear J1 J2

        % action_range1 = -3:0.1:3;
        % 
        % plot_cost_function(task1,action_range);
        % plot_cost_function(task2,action_range);
        % 
        convergence(game)
        %% Define the model
        % SENSORY SYSTEM
        % sensory input equals partner action + own action + noise

        H1 = [eye(xsz1)];
        H2 = [eye(xsz2)];

        L1 = [-0.*eye(xsz1)];
        L2 = [-0.*eye(xsz2)];

        % sensory noise
        sigmay1 = sigmay_other(g)*eye(xsz1);
        sigmay2 = sigmay_other(g)*eye(xsz2);

        % PARTNER MODEL
        % Retention rate
        A1 = A(g)*eye(xsz1);
        A2 = A(g)*eye(xsz1);

        % Process oise
        sigmax1 = sigmax(g)*eye(xsz1); %
        sigmax2 = sigmax(g)*eye(xsz2); %

        % Initial state mean and variance
        x01 = X0(g).*ones(xsz1,1);
        P01 = P0(g).*eye(xsz1);
        x02 = X0(g).*ones(xsz2,1);
        P02 = P0(g).*eye(xsz2);

        % ACTION CONTROLLER
        % Temperature d0ecreases during experiment
        lambda1_ini = lambda(g);%.2;
        decay_rate1 = lambda_decay(g);
        lambda2_ini = lambda(g);
        decay_rate2 = lambda_decay(g);

        % INITIALIZATION
        % Initialize the observer(s)
        observer1 = observer(A1,[],H1,L1,sigmax1,sigmay1,x01,P01); % V2 considering also own action in sensory input
        % observer1 = observer(A1,[],H1,[],sigmax1,sigmay1,x01,P01);
        %K1stat = StationaryK(observer1);

        observer2 = observer(A2,[],H2,L2,sigmax2,sigmay2,x02,P02);
        % observer2 = observer(A2,[],H2,[],sigmax2,sigmay2,x02,P02);
        %K2stat = StationaryK(observer2);

        % Initialize the simulation
        P1_prior = cell(T,1);
        P2_prior = cell(T,1);
        x1_prior = zeros(xsz1,T);
        x2_prior = zeros(xsz2,T);

        x1_prior(:,1) = x01;
        P1_prior{1} = P01;

        x2_prior(:,1) = x02;
        P2_prior{1} = P02;

        u1final = [];
        u2final = [];

        lambda1(1) = lambda1_ini; % comment to start over the controller at each block
        lambda2(1) = lambda2_ini;
        clear contr1 contr2
        contr1 = controller(task1,lambda1_ini,decay_rate1,action_range(1),action_range(end));
        contr2 = controller(task2,lambda2_ini,decay_rate2,action_range(1),action_range(end));
        close all

        for b = 1:n_blocks

            clear pu1_x1 pu2_x2

            block_trials = ((b-1)*n_reps+1):b*n_reps;

            for t= block_trials
                t
                
                % Generate new actions
                [u1(:,t),decision1(:,t), contr1] = generate_action(contr1, x1_prior(:,t),P1_prior{t},lambda1(t));
                [u2(:,t),decision2(:,t), contr2] = generate_action(contr2, x2_prior(:,t),P2_prior{t},lambda2(t));


                if check_single_steps
                    figure('pos',[0 0 200 100])
                    xline(u1(:,t),'LineWidth',2,'Color','k')
                    xlim([action_range(1) action_range(end)])
                    yticks([])
                    xticks([])

                    figure('pos',[0 100 200 100])
                    xline(u2(:,t),'LineWidth',2,'Color','k')
                    xlim([action_range(1) action_range(end)])
                    yticks([])
                    xticks([])
                end

                jj1_temp = cost(task1,u1(:,t),x1_prior(:,t));

                jj1r = jj1_temp{1};
                jj1s = jj1_temp{2};

                jj2_temp = cost(task2,u2(:,t),x2_prior(:,t));
                jj2r = jj2_temp{1};
                jj2s = jj2_temp{2};

                pu1r_x1(:,t) = exp(-jj1r./lambda1(t));
                pu2r_x2(:,t) = exp(-jj2r./lambda2(t));

                pu1s_x1(:,t) = exp(-jj1s./lambda1(t));
                pu2s_x2(:,t) = exp(-jj2s./lambda2(t));

                prior1(t,:) = get(contr1,'prior');
                prior2(t,:) = get(contr2,'prior');

                % Generate sensory feedback
                % y1(:,t) = generate_sensory(observer1,u2imposed(:,t),u1(:,t),'sim');
                % y2(:,t) = generate_sensory(observer2,u1imposed(:,t),u2(:,t),'sim');

                y1(:,t) = generate_sensory(observer1,u2(:,t),u1(:,t),'sim');
                y2(:,t) = generate_sensory(observer2,u1(:,t),u2(:,t),'sim');

                if check_single_steps
                    y_norm = normpdf(action_range,H1*u2(:,t) + L1*u1(:,t),sigmay1)
                     y_norm = normpdf(action_range,H1*u2(:,t),sigmay1)
                    figure('pos',[0 200 200 100])
                    plot(action_range,y_norm,'LineWidth',2,'Color','r')
                    xline(y1(:,t),'LineWidth',2,'Color','r')
                    xlim([action_range(1) action_range(end)])
                    yticks([])
                    box off
                    xticks([])
                end

                % Update state observers: static
                %[x1_post(:,t),P1_post{t},x1_prior(:,t+1),P1_prior{t+1}] = kalman_estimate_stationary(observer1,x1_prior(:,t),P1_prior{t},y1(:,t),K1stat);
                %[x2_post(:,t),P2_post{t},x2_prior(:,t+1),P2_prior{t+1}] = kalman_estimate_stationary(observer2,x2_prior(:,t),P2_prior{t},y2(:,t),K2stat);

                % Update state observers: dynamic
                [x1_post(:,t),P1_post{t},x1_prior(:,t+1),P1_prior{t+1},Kalman1{t}] = kalman_estimate(observer1,x1_prior(:,t),P1_prior{t},y1(:,t),u1(:,t));
                [x2_post(:,t),P2_post{t},x2_prior(:,t+1),P2_prior{t+1},Kalman2{t}] = kalman_estimate(observer2,x2_prior(:,t),P2_prior{t},y2(:,t),u2(:,t));
                %     [x1_post(:,t),P1_post{t},x1_prior(:,t+1),P1_prior{t+1},Kalman1{t}] = kalman_estimate(observer1,x1_prior(:,t),P1_prior{t},y1(:,t),[]);
                %     [x2_post(:,t),P2_post{t},x2_prior(:,t+1),P2_prior{t+1},Kalman2{t}] = kalman_estimate(observer2,x2_prior(:,t),P2_prior{t},y2(:,t),[]);

                % Update temperature
                lambda1(t+1) = decay_rate1*lambda1(t);
                lambda2(t+1) = decay_rate2*lambda2(t);

                 if check_single_steps
                    x_norm = normpdf(action_range,x1_prior(:,t+1),sigmax1)
                    
                    figure('pos',[0 300 200 100])
                    plot(action_range,x_norm,'LineWidth',2,'Color','b')
                    % xline(x1_prior(:,t),'--','LineWidth',2,'Color','b')
                    xlim([action_range(1) action_range(end)])
                    yticks([])
                    box off
                    xticks([])

                    for i = 1:length(action_range)
                        jj1_pred = cost(task1,action_range(i),x1_prior(:,t+1));
                        jj1_pred_r(i) = jj1_pred{1};
                        jj1_pred_s(i) = jj1_pred{2};
                    end

                    figure('pos',[0 400 200 100])
                    plot(action_range,jj1_pred_r,'LineWidth',2,'Color',[174 5 54]./255)
                    hold on
                    plot(action_range,jj1_pred_s,'LineWidth',2,'Color',[131 4 41]./255)
                    box off
                    yticks([])
                    xlim([action_range(1) action_range(end)])

                    [mr,mridx] = min(jj1_pred_r);
                    [ms,msidx] = min(jj1_pred_s);

                    plot(action_range(mridx),jj1_pred_r(mridx),'*','Color',[174 5 54]./255,'linewidth',2)
                    plot(action_range(msidx),jj1_pred_s(msidx),'*','Color',[131 4 41]./255,'linewidth',2)
                    xticks([])


                    if mr < ms
                        u_norm = normpdf(action_range,action_range(mridx),lambda1(t+1))
                    else
                        u_norm = normpdf(action_range,action_range(msidx),lambda1(t+1))
                    end

                    figure('pos',[0 500 200 100])
                    plot(action_range,u_norm,'LineWidth',2,'Color','k')
                    xlim([action_range(1) action_range(end)])
                    yticks([])
                    xticks([])
                    box off
                    
                 end

            end
            final_trials = (b-1)*n_reps+((round(n_reps/2)+1):n_reps);%(b-1)*n_reps+(21:40);

            u1final = [u1final u1(:,final_trials)];
            u2final = [u2final u2(:,final_trials)];
        end
        %% Plot Decisions and Actions

        plot_action_decision(sh_game,u1,u2,decision1,decision2,action_range)
        savefig(fullfile(figdir,['actions_decisions' num2str(n_sim) '.fig']));

        plot_action(sh_game,u1,u2,1:T,action_range)
        ylim([-1.5 1.5])
        xlim([-1.5 1.5])

        plot_action(sh_game,u1final,u2final,1:T,action_range)
        ylim([-1.5 1.5])
        xlim([-1.5 1.5])

        savefig(fullfile(figdir,['actions' num2str(n_sim) '.fig']));

        %% Plot lambdas
        % figure
        % set(gcf,'pos',[100 100 300 150])
        % hold on
        % plot(1:T,lambda1(1:T),'b')
        % plot(1:T,lambda2(1:T),'r')
        % xlabel('Trials')
        % ylabel('\lambda_i')

        %% Probabilities of NEs and non-cooperative solution

        for t = 1:T
            SS(d,t) = u1(t)>mean(action_range) && u2(t)>mean(action_range);
            HH(d,t) = u1(t)<mean(action_range) && u2(t)<mean(action_range);
            SH(d,t) = u1(t)>mean(action_range) && u2(t)<mean(action_range);
            HS(d,t) = u1(t)<mean(action_range) && u2(t)>mean(action_range);
            
        end

        for b = 1:n_blocks
            curr_reps = (b-1)*n_reps + (1:n_reps);
            pSS{g}(d,b) = mean(SS(d,curr_reps)); seSS{g}(d,b) = std(SS(d,curr_reps))./sqrt(n_reps);
            pHH{g}(d,b) = mean(HH(d,curr_reps)); seHH{g}(d,b) = std(HH(d,curr_reps))./sqrt(n_reps);
            pSH{g}(d,b) = mean(SH(d,curr_reps)); seSH{g}(d,b) = std(SH(d,curr_reps))./sqrt(n_reps);
            pHS{g}(d,b) = mean(HS(d,curr_reps)); seHS{g}(d,b) = std(HS(d,curr_reps))./sqrt(n_reps);
        end

        % probabilità condizionate
        for t = 2:T

            S1_S2(t) = u1(t)>mean(action_range) && u2(t-1)>mean(action_range);
            S2_S1(t) = u2(t)>mean(action_range) && u1(t-1)>mean(action_range);
            S1_H2(t) = u1(t)>mean(action_range) && u2(t-1)<mean(action_range);
            S2_H1(t) = u2(t)>mean(action_range) && u1(t-1)<mean(action_range);

            H1_S2(t) = u1(t)<mean(action_range) && u2(t-1)>mean(action_range);
            H2_S1(t) = u2(t)<mean(action_range) && u1(t-1)>mean(action_range);
            H1_H2(t) = u1(t)<mean(action_range) && u2(t-1)<mean(action_range);
            H2_H1(t) = u2(t)<mean(action_range) && u1(t-1)<mean(action_range);

        end

        for b = 1:n_blocks
            curr_reps = (b-1)*n_reps + (1:n_reps);
            pS1_S2{g}(d,b) = mean(S1_S2(curr_reps)); seS1_S2{g}(d,b) = std(S1_S2(curr_reps))./sqrt(n_reps);
            pS2_S1{g}(d,b) = mean(S2_S1(curr_reps)); seS2_S1{g}(d,b) = std(S2_S1(curr_reps))./sqrt(n_reps);
            pS1_H2{g}(d,b) = mean(S1_H2(curr_reps)); seS1_S2{g}(d,b) = std(S1_H2(curr_reps))./sqrt(n_reps);
            pS2_H1{g}(d,b) = mean(S2_H1(curr_reps)); seS2_S1{g}(d,b) = std(S2_H1(curr_reps))./sqrt(n_reps);

            pH1_S2{g}(d,b) = mean(H1_S2(curr_reps)); seH1_S2{g}(d,b) = std(H1_S2(curr_reps))./sqrt(n_reps);
            pH2_S1{g}(d,b) = mean(H2_S1(curr_reps)); seH2_S1{g}(d,b) = std(H2_S1(curr_reps))./sqrt(n_reps);
            pH1_H2{g}(d,b) = mean(H1_H2(curr_reps)); seH1_S2{g}(d,b) = std(H1_H2(curr_reps))./sqrt(n_reps);
            pH2_H1{g}(d,b) = mean(H2_H1(curr_reps)); seH2_S1{g}(d,b) = std(H2_H1(curr_reps))./sqrt(n_reps);
        end

        %if d == 1
            save(fullfile(cd,'results',['sim' num2str(n_sim) '_D' num2str(d)]),'u1','u2','u1imposed','u2imposed',...
                'A1','A2','sigmay1','sigmay2','sigmax1','sigmax2','x01','x02','P01','P02','lambda1_ini','lambda2_ini','decay_rate1','decay_rate2')
        %end
    end
    
    %%
    for s = 1:n_sub
        figure
        subplot(211)
        hold on
        patch([1:n_blocks n_blocks:-1:1],[pSS{g}(s,:)-seSS{g}(s,:) fliplr(pSS{g}(s,:)+seSS{g}(s,:))],purple,'FaceAlpha',0.1,'EdgeAlpha',0)
        patch([1:n_blocks n_blocks:-1:1],[pHH{g}(s,:)-seHH{g}(s,:) fliplr(pHH{g}(s,:)+seHH{g}(s,:))],orange,'FaceAlpha',0.1,'EdgeAlpha',0)
        c(1) = plot(pSS{g}(s,:),'-','color',purple,'linewidth',2)
        c(2) = plot(pHH{g}(s,:),'-','color',orange,'linewidth',2)
        legend(c,{'SS','RR'},'fontname','times')
        box off
        xlabel('Blocks','fontname','times')
        ylabel('Probability','fontname','times')
        legend boxoff
        ylim([0 1])

        subplot(212)
        hold on
        patch([1:n_blocks n_blocks:-1:1],[pSH{g}(s,:)-seSH{g}(s,:) fliplr(pSH{g}(s,:)+seSH{g}(s,:))],light_green,'FaceAlpha',0.1,'EdgeAlpha',0)
        patch([1:n_blocks n_blocks:-1:1],[pHS{g}(s,:)-seHS{g}(s,:) fliplr(pHS{g}(s,:)+seHS{g}(s,:))],green,'FaceAlpha',0.1,'EdgeAlpha',0)
        c(1) = plot(pSH{g}(s,:),'-','color',light_green,'linewidth',2)
        c(2) = plot(pHS{g}(s,:),'-','color',green,'linewidth',2)
        legend(c,{'SR','RS'},'fontname','times')
        box off
        legend boxoff
        xlabel('Blocks','fontname','times')
        ylabel('Probability','fontname','times')
        ylim([0 1])
        set(gcf,'Position',[0 300 150 300])
        savefig(fullfile(figdir,['prob_actions' num2str(n_sim) 'd' num2str(s) '.fig']));
    end
    
    %% pop level P(HH),P(SS)...
    
    if n_sub > 1
        figure
        subplot(211)
        hold on
        patch([1:n_blocks n_blocks:-1:1],[mean(pSS{g})-std(pSS{g})./sqrt(n_sub) fliplr(mean(pSS{g})+std(pSS{g})./sqrt(n_sub))],purple,'FaceAlpha',0.1,'EdgeAlpha',0)
        patch([1:n_blocks n_blocks:-1:1],[mean(pHH{g})-std(pHH{g})./sqrt(n_sub) fliplr(mean(pHH{g})+std(pHH{g})./sqrt(n_sub))],orange,'FaceAlpha',0.1,'EdgeAlpha',0)
        c(1) = plot(mean(pSS{g}),'-','color',purple,'linewidth',2)
        c(2) = plot(mean(pHH{g}),'-','color',orange,'linewidth',2)
        yline(0.25)
        legend(c,{'SS','RR'},'fontname','times')
        box off
        xlabel('Blocks','fontname','times')
        ylabel('Probability','fontname','times')
        set(gca,'FontName','times')
        xlim([0 n_blocks+1])
        legend boxoff
        ylim([0 1])

        subplot(212)
        hold on
        patch([1:n_blocks n_blocks:-1:1],[mean(pSH{g})-std(pSH{g})./sqrt(n_sub) fliplr(mean(pSH{g})+std(pSH{g})./sqrt(n_sub))],light_green,'FaceAlpha',0.1,'EdgeAlpha',0)
        patch([1:n_blocks n_blocks:-1:1],[mean(pHS{g})-std(pHS{g})./sqrt(n_sub) fliplr(mean(pHS{g})+std(pHS{g})./sqrt(n_sub))],green,'FaceAlpha',0.1,'EdgeAlpha',0)
        c(1) = plot(mean(pSH{g}),'-','color',light_green,'linewidth',2)
        c(2) = plot(mean(pHS{g}),'-','color',green,'linewidth',2)
        yline(0.25)
        legend(c,{'SR','RS'},'fontname','times')
        box off
        xlim([0 n_blocks+1])
        legend boxoff
        xlabel('Blocks','fontname','times')
        ylabel('Probability','fontname','times')
        ylim([0 1])
        set(gca,'FontName','times')
        
        set(gcf,'Position',[0 300 150 300])
        savefig(fullfile(figdir,['prob_actions' num2str(n_sim) 'd' num2str(s) '.fig']));
    end
    %% pop level P(H1tH2t-1),P(S1tS2t-1)...
    if n_sub > 1
        figure
        subplot(221)
        hold on
        patch([1:n_blocks n_blocks:-1:1],[mean(pS1_S2{g})-std(pS1_S2{g})./sqrt(n_sub) fliplr(mean(pS1_S2{g})+std(pS1_S2{g})./sqrt(n_sub))],purple,'FaceAlpha',0.1,'EdgeAlpha',0)
        patch([1:n_blocks n_blocks:-1:1],[mean(pS1_H2{g})-std(pS1_H2{g})./sqrt(n_sub) fliplr(mean(pS1_H2{g})+std(pS1_H2{g})./sqrt(n_sub))],orange,'FaceAlpha',0.1,'EdgeAlpha',0)
        c(1) = plot(mean(pS1_S2{g}),'-','color',purple,'linewidth',2)
        c(2) = plot(mean(pS1_H2{g}),'-','color',orange,'linewidth',2)
        legend(c,{'p(S1|S2)','p(S1|R2)'},'fontname','times')
        legend box off

        subplot(222)
        hold on
        patch([1:n_blocks n_blocks:-1:1],[mean(pS2_S1{g})-std(pS2_S1{g})./sqrt(n_sub) fliplr(mean(pS2_S1{g})+std(pS2_S1{g})./sqrt(n_sub))],purple,'FaceAlpha',0.1,'EdgeAlpha',0)
        patch([1:n_blocks n_blocks:-1:1],[mean(pS2_H1{g})-std(pS2_H1{g})./sqrt(n_sub) fliplr(mean(pS2_H1{g})+std(pS2_H1{g})./sqrt(n_sub))],orange,'FaceAlpha',0.1,'EdgeAlpha',0)
        c(1) = plot(mean(pS2_S1{g}),'-','color',purple,'linewidth',2)
        c(2) = plot(mean(pS2_H1{g}),'-','color',orange,'linewidth',2)
        legend(c,{'p(S2|S1)','p(S2|R2)'},'fontname','times')
        box off
        xlabel('Blocks','fontname','times')
        ylabel('Probability','fontname','times')
        legend boxoff
        ylim([0 1])

        subplot(223)
        hold on
        patch([1:n_blocks n_blocks:-1:1],[mean(pH1_S2{g})-std(pH1_S2{g})./sqrt(n_sub) fliplr(mean(pH1_S2{g})+std(pH1_S2{g})./sqrt(n_sub))],purple,'FaceAlpha',0.1,'EdgeAlpha',0)
        patch([1:n_blocks n_blocks:-1:1],[mean(pH1_H2{g})-std(pH1_H2{g})./sqrt(n_sub) fliplr(mean(pH1_H2{g})+std(pH1_H2{g})./sqrt(n_sub))],orange,'FaceAlpha',0.1,'EdgeAlpha',0)
        c(1) = plot(mean(pH1_S2{g}),'-','color',purple,'linewidth',2)
        c(2) = plot(mean(pH1_H2{g}),'-','color',orange,'linewidth',2)
        legend(c,{'p(R1|S2)','p(R1|R2)'},'fontname','times')
        box off
        legend boxoff
        xlabel('Blocks','fontname','times')
        ylabel('Probability','fontname','times')
        ylim([0 1])

        subplot(224)
        hold on
        patch([1:n_blocks n_blocks:-1:1],[mean(pH2_S1{g})-std(pH2_S1{g})./sqrt(n_sub) fliplr(mean(pH2_S1{g})+std(pH2_S1{g})./sqrt(n_sub))],purple,'FaceAlpha',0.1,'EdgeAlpha',0)
        patch([1:n_blocks n_blocks:-1:1],[mean(pH2_H1{g})-std(pH2_H1{g})./sqrt(n_sub) fliplr(mean(pH2_H1{g})+std(pH2_H1{g})./sqrt(n_sub))],orange,'FaceAlpha',0.1,'EdgeAlpha',0)
        c(1) = plot(mean(pH2_S1{g}),'-','color',purple,'linewidth',2)
        c(2) = plot(mean(pH2_H1{g}),'-','color',orange,'linewidth',2)
        legend(c,{'p(R2|S1)','p(R2|R1)'},'fontname','times')
        box off
        legend boxoff
        xlabel('Blocks','fontname','times')
        ylabel('Probability','fontname','times')
        ylim([0 1])

        set(gcf,'Position',[0 300 300 300])
        savefig(fullfile(figdir,['prob_cond_actions' num2str(n_sim) 'd' num2str(s) '.fig']));
    end
end

save(fullfile(cd,'results',['sim' num2str(n_sim)]),"sigmax",'sigmay_other','A','lambda','lambda_decay','P0','X0',...
    'SS','HH','SH','HS','pSS','pHH','pSH','pHS','seSS','seHH','seSH','seHS')
