% Fit SH Game
close all
clear all
clear classes

addpath("..");

datadir = 'results';
n_sub = 8;
for d = 1:n_sub
    fname = ['sim12_D' num2str(d) '.mat']; % d12 is the simulation with variability

    load(fullfile(datadir,fname));

    %% DEFINE THE TASK
    action_range = [-1:0.01:1];%[0:0.01:1];%
    uR = -1; % Rabbit location
    uS = 1; % Stag location
    u01 = 0;%mean(action_range); % Initial position 1
    u02 = 0;%mean(action_range); %Initial position 2

    % % Gain term
    w = 9/2;
    wR = w; wS = w;
    zS = 1;%4/9;
    zR = 5;%20/9;

    sh_game = stagHunt_game(uR, uS, zR, zS, u01,u02,wR,wS);
    task1 = get(sh_game,'task1')
    task2 = get(sh_game,'task2')

    xsz1 = get(task1,'xsize');
    xsz2 = get(task2,'xsize');

    % DISPLAY Nash equilibria
    game = quadraticgame(task1,task2);

    U1 = u1';
    U1imp = u1';%u1imposed';%
    U2 = u2';
    U2imp = u2';%u2imposed';%

    %% Define player
    % Initial values of the parameters

    parnames = {'A','B','C','D','SigmaX','SigmaY','x0','P0','lambda0','a'};

    % X0 = [(-1/log(0.9)),  .1, .1,      eps,   2*eps,   .1,  (-1/log(0.9))];
    % LB = [(-1/log(0.8)),  eps, eps, -1e-10,     0,  0.09,  (-1/log(0.8))];
    % UB = [(-1/log(0.9999)), 10.1, 10.1,  1e-10, 1e-10,  0.11, (-1/log(0.9999))];

    X0 = [0.9,  5, 5,      eps,   2*eps,   .1,  0.9];
    LB = [0.8,  eps, eps, -1e-10,     0,  0.09,  0.8];
    %UB = [0.9999, 10.1, 10.1,  1e-10, 1e-10,  0.11, 0.9999];
    UB = [0.9999, 11, 11,  1e-10, 1e-10,  0.11, 0.9999];


    opts = optimoptions('fmincon','FunctionTolerance',1e-15);
    opts.Display = 'iter';

    f = @(x) costfunction(x,task1,U2,U1);
    f(X0)

    Aeq = [0,0,0,1,0,0,0];
    beq = [1];

    % partner 1
    X1(d,:) = fmincon(@(x) costfunction(x,task1,U2imp,U1),X0,[],[],[],[],LB,UB,[],opts);

    % partner 2
    X2(d,:) = fmincon(@(x) costfunction(x,task2,U1imp,U2),X0,[],[],[],[],LB,UB,[],opts);

    [X0' X1(d,:)' X2(d,:)']



    % %%
    % X1True =    [-1/log(A1) sigmax1 sigmay1 x01 P01 lambda1_ini -1/log(decay_rate1)];
    % pars1True = {exp(-1/X1True(1))*eye(xsz1),[],eye(xsz1),eye(xsz1),X1True(2)*eye(xsz1),X1True(3)*eye(xsz1),X1True(3+(1:xsz1)),X1True(3+xsz1+1)*eye(xsz1),X1True(3+xsz1+2),exp(-1/X1True(3+xsz1+3))};
    % model1True = playermodel(task1,pars1True,parnames);
    % %[cost1True,U1hatTrue,Phat_U1True,pi1True] = pseudolik(model1True,U2imp',U1');
    % cost1True =costfunction(X1True,task1,U2imp,U1);
    %
    % X2True = [-1/log(A2) sigmax2 sigmay2 x01 P02 lambda2_ini -1/log(decay_rate2)];
    % pars2True = {exp(-1/X2True(1))*eye(xsz2),[],eye(xsz2),eye(xsz2),X2True(2)*eye(xsz2),X2True(3)*eye(xsz2),X2True(3+(1:xsz2)),X2True(3+xsz1+1)*eye(xsz2),X2True(3+xsz2+2),exp(-1/X2True(3+xsz2+3))};
    % model2True = playermodel(task2,pars2True,parnames);
    % %[cost2True,U2hatTrue,Phat_U2True,pi2True] = pseudolik(model2True,U1imp',U2');
    % cost2True =costfunction(X2True,task2,U1imp,U2);
    %
    % parnames = {'A','B','C','D','SigmaX','SigmaY','x0','P0','lambda0','a'};
    % pars1 = {exp(-1/X1(1))*eye(xsz1),[],eye(xsz1),eye(xsz1),X1(2)*eye(xsz1),X1(3)*eye(xsz1),X1(3+(1:xsz1)),X1(3+xsz1+1)*eye(xsz1),X1(3+xsz1+2),exp(-1/X1(3+xsz1+3))};
    % model1 = playermodel(task1,pars1,parnames);
    % [cost1,U1hat,Phat_U1,pi1] = pseudolik(model1,U2imp',U1');
    % cost1 = costfunction(X1,task1,U2imp,U1);
    %
    % pars2 = {exp(-1/X2(1))*eye(xsz2),[],eye(xsz2),eye(xsz2),X2(2)*eye(xsz2),X2(3)*eye(xsz2),X2(3+(1:xsz2)),X2(3+xsz1+1)*eye(xsz2),X2(3+xsz2+2),exp(-1/X2(3+xsz2+3))};
    % model2 = playermodel(task2,pars2,parnames);
    % [cost2,U2hat,Phat_U2,pi2] = pseudolik(model2,U1imp',U2');
    % cost2 =costfunction(X2,task2,U1imp,U2);
    %
    % cost1True>cost1
    %
    % cost2True>cost2

    %%
    X1True(d,:) =    [A1 sigmax1 sigmay1 x01 P01 lambda1_ini decay_rate1];
    pars1True = {X1True(d,1)*eye(xsz1),[],eye(xsz1),eye(xsz1),X1True(d,2)*eye(xsz1),X1True(d,3)*eye(xsz1),X1True(d,3+(1:xsz1)),X1True(d,3+xsz1+1)*eye(xsz1),X1True(d,3+xsz1+2),X1True(d,3+xsz1+3)};
    model1True = playermodel(task1,pars1True,parnames);
    %[cost1True,U1hatTrue,Phat_U1True,pi1True] = pseudolik(model1True,U2imp',U1');
    cost1True =costfunction(X1True(d,:),task1,U2imp,U1);

    X2True(d,:) = [A2 sigmax2 sigmay2 x01 P02 lambda2_ini decay_rate2];
    pars2True = {X2True(d,1)*eye(xsz2),[],eye(xsz2),eye(xsz2),X2True(d,2)*eye(xsz2),X2True(d,3)*eye(xsz2),X2True(d,3+(1:xsz2)),X2True(d,3+xsz1+1)*eye(xsz2),X2True(d,3+xsz2+2),X2True(d,3+xsz2+3)};
    model2True = playermodel(task2,pars2True,parnames);
    %[cost2True,U2hatTrue,Phat_U2True,pi2True] = pseudolik(model2True,U1imp',U2');
    cost2True =costfunction(X2True(d,:),task2,U1imp,U2);

    parnames = {'A','B','C','D','SigmaX','SigmaY','x0','P0','lambda0','a'};
    pars1 = {X1(d,1)*eye(xsz1),[],eye(xsz1),eye(xsz1),X1(d,2)*eye(xsz1),X1(d,3)*eye(xsz1),X1(d,3+(1:xsz1)),X1(d,3+xsz1+1)*eye(xsz1),X1(d,3+xsz1+2),X1(d,3+xsz1+3)};
    model1 = playermodel(task1,pars1,parnames);
    [cost1,U1hat,Phat_U1,pi1] = pseudolik(model1,U2imp',U1');
    cost1 = costfunction(X1(d,:),task1,U2imp,U1);

    pars2 = {X2(d,1)*eye(xsz2),[],eye(xsz2),eye(xsz2),X2(d,2)*eye(xsz2),X2(d,3)*eye(xsz2),X2(d,3+(1:xsz2)),X2(d,3+xsz1+1)*eye(xsz2),X2(d,3+xsz2+2),X2(d,3+xsz2+3)};
    model2 = playermodel(task2,pars2,parnames);
    [cost2,U2hat,Phat_U2,pi2] = pseudolik(model2,U1imp',U2');
    cost2 = costfunction(X2(d,:),task2,U1imp,U2);

    cost1True>cost1
    cost2True>cost2

    %% Plotting weighted action
    for t = 1:size(Phat_U1,1)
        uhat1_t(:,t) = zeros(1,1);
        uhat2_t(:,t) = zeros(1,1);

        for s = 1:size(pi1,2)
            uhat1_t(:,t) = uhat1_t(:,t) + pi1(t,s).*U1hat{t,s};
            uhat2_t(:,t) = uhat2_t(:,t) + pi2(t,s).*U2hat{t,s};
        end
    end

    %% Correlation between true actions and estimated action
    c1temp = corrcoef(U1,uhat1_t');
    c1(d) = c1temp(1,2).^2;
    c2temp = corrcoef(U2,uhat2_t');
    c2(d) = c2temp(1,2).^2;

    %% Action Visualization
    figure
    set(gcf,'pos',[0 200 500 200])
    subplot(211)
    plot(U1,'-k','linewidth',1.5)
    hold on
    plot(uhat1_t,':b','linewidth',1.5)
    title(['$\hat{u_1}$ $R^2$ = ' num2str(c1)],'interpreter','latex')
    %title(['Estimated Actions Player1 $R^2$ = ' num2str(c1)]);%,'interpreter','latex')
    yticks([-1 0 1])
    ylabel('Actions','interpreter','latex')
    %yticklabels({'R','S'})
    box off
    set(gca,'fontname','times')


    subplot(212)
    plot(U2,'-k','linewidth',1.5)
    hold on
    plot(uhat2_t,':r','linewidth',1.5)
    title(['$\hat{u_2}$ $R^2$ = ' num2str(c2)],'interpreter','latex')
    %title(['Estimated Actions Player2 $R^2$ = ' num2str(c2)]);%,'interpreter','latex')
    xlabel('Trials','interpreter','latex')
    ylabel('Actions','interpreter','latex')
    yticks([-1 0 1])
    %yticklabels({'R','S'})
    box off
    set(gca,'fontname','times')
    set(gcf,'pos',[300 0 180 250])
    set(gcf,'pos',[300 0 350 160])


end
%% Parameters percentage variation
AvErrParams1 = mean((X1-X1True));
AvErrParams2 = mean((X2-X2True));
AvErrParams = mean([X1; X2]-[X1True; X2True]);

SErrParams1 = std((X1-X1True))./sqrt(n_sub);
SErrParams2 = std((X2-X2True))./sqrt(n_sub);
varnames = {'A','SigmaX','SigmaY','x0','P0','lambda0','a'};

Pl1_Table = table(mean(X1)', mean(X1True)', AvErrParams1');%,'RowNames',varnames);
Pl1_Table.Properties.VariableNames = ["Average Estimated","True","Average Error"];
Pl2_Table = table(mean(X2)', mean(X2True)', AvErrParams2');
Pl2_Table.Properties.VariableNames = ["Average Estimated","True","Average Error"];
Pl1_Table
Pl2_Table
av_percVarP1 = mean(AvErrParams1);
av_percVarP2 = mean(AvErrParams2);

Pl_Table = table(mean([X1;X2])', mean([X1True;X2True])', AvErrParams');
Pl_Table.Properties.VariableNames = ["Average Estimated Parameters","True Parameters","Average Error"];
Pl_Table.Properties.RowNames = varnames'

mean([c1 c2])
