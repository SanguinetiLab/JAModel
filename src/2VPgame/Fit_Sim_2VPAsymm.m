close all
clear all
clc

v = input('Set version number: ');

suffix = ['fit2VPAsymm_multistrat_V' num2str(v)]; %

% Flags
fit = 0; % 0 if you wanto just to check parameters and simulate
check = 0;
simulate_single_dyads = 1;
simulate_global_groups = 0;
mode = 'dyad';% 'single';
gen_action_mode = 'stochastic';%'deterministic';%
n_sim = 1; % number of simulations per dyad

%% Getting Dataset
cd ../..
rootfolder = cd;
datadir = fullfile(rootfolder,'data','experiments','expAsymVPs');
cd src/2VPgame/  

% Set the numebr of nodes
nnodes_dataset = 5; %10

%% Exp Design
n_groups = 3;
n_dyads_group = 5;
n_dyads = n_groups*n_dyads_group;
subH = 1:5; subVH = 6:10; subPV = 11:15;

%% Loading forces to remove unconnected trials, tc11 and tc22, ts1, ts2, te1, te2
load(fullfile(datadir, ['forces.mat']));

%% Fitting for each Dyad
if fit
    for d = 1:n_dyads
        % Reading
        if nnodes_dataset== 5
            data_action = readtable(fullfile(datadir,['dyad' num2str(d) '_5nodes.dat']));
        else
            data_action = readtable(fullfile(datadir,['dyad' num2str(d) '.dat']));
        end

        % Removing unconnected trials
        unconnected = find(forces(d,:) == 0);
        data_action(unconnected,:) = [];
        nodes_num = size(data_action,2)/4; % number of nodes in spline approximation

        % Partner 1 actions
        U1 = table2array(data_action(:,1:(2*nodes_num)));
        % Partner 2 actions
        U2 = table2array(data_action(:,2*nodes_num+(1:(2*nodes_num))));

        if nnodes_dataset == 5
            U1(:,1:2:9) = U1(:,1:2:9)-0.05;
            U2(:,1:2:9) = U2(:,1:2:9)-0.05;
        end

        % Defining task specs
        [tcE,tcL, duration, yrange, start, final, VP1, VP2, tgterror,VPerror] = get_2VPtask_specs('Asymm');

        T = size(data_action,1); % number of plays

        % Define 2-VP game
        twoVPg = twoVP_game(nodes_num,VP1,VP2,tcE,tcL,start,final,duration,VPerror,tgterror);
        task1 = get(twoVPg,'task1');
        task2 = get(twoVPg,'task2');
        xsz1 = get(task1,'xsize');
        xsz2 = get(task2,'xsize');

        game = quadraticgame(task1,task2);

        % Define Players
        % Initial values of the parameters

        % Mu0 initialization
        x0 = zeros(xsz1,1);%)eps.*ones(xsz1,1);
        dt = duration/(nodes_num+2);
        tt = 0:dt:duration
        tt_half = 0:dt:(duration/2);
        for n = 1:nodes_num
            xmin_jeark(n) =  start(1) + (final(1)-start(1))*minimum_jerk_traj(tt(n+1)/duration);
            x0(n*2-1) =  xmin_jeark(n) ;
            x0_LB((n*2-1):n*2) = [xmin_jeark(n)-0.0001 -0.0001];%[xmin_jeark(n)-0.01 -0.01];%yrange(1)/2];%[0.1*x01(n*2-1) -0.05];
            x0_UB((n*2-1):n*2) = [xmin_jeark(n)+0.0001 0.0001];%[xmin_jeark(n)+0.01 0.01];%yrange(2)/2];%[10*x01(n*2-1) 0.05];
        end

        % V15
        Xtrue = [0.9,  0.001,  0.01,       x0',      0.0001,  0.5,    0.99];
        LB =    [0.85,  eps,    eps,    x0_LB,    eps,     eps,    0.9];
        UB =    [1-eps, 0.03,   0.06,      x0_UB,    0.03,   1,      1-eps];
        if d == 7 || d == 8 || d == 10 || d == 13
            Xtrue = [0.999,  0.001,  0.01,       x0',      0.0001,  0.5,    0.999];
        end

        % V16
        Xtrue = [0.9,  0.001,  0.01,       x0',      0.0001,  0.5,    0.99];
        LB =    [eps,  eps,    eps,    x0_LB,    eps,     eps,   eps];
        UB =    [1-eps, 0.03,   0.06,      x0_UB,    0.03,   10,      1-eps];

        % V18
        Xtrue = [0.9,  0.001,  0.01,       x0',      0.0001,  0.5,    0.99];
        LB =    [eps,  eps,    eps,    x0_LB,    eps,     eps,   eps];
        UB =    [0.9999, 1e-5,   0.06,      x0_UB,    1e-5,   5,     0.9999];
        if d == 2 || d ==7 
            Xtrue = [0.99,  0.001,  0.01,       x0',      0.0001,  0.5,    0.99];
        elseif d == 3 || d == 5 
            Xtrue = [0.999,  0.001,  0.01,       x0',      0.0001,  0.5,    0.999];
        elseif d == 6 || d == 8 || d ==9 || d ==10 || d ==11 || d ==12 || d == 13 || d == 14
            Xtrue = [0.999,  0.001,  0.01,       x0',      0.00001,  1,    0.999];
        end

        % V19
        Xtrue = [0.9,  1e-4,  0.01,       x0',      1e-8,  eps,    0.99];
        LB =    [eps,  eps,    eps,    x0_LB,    eps,     eps,   eps];
        UB =    [0.999, 1e-4,   0.2,      x0_UB,    1e-7,   10,     0.999];
        if d == 8 || d == 13
            Xtrue = [0.99,  1e-4,  0.01,       x0',      1e-8,  eps,    0.99];

        end

        % V20 Like V7
        Xtrue = [0.99,0.05^2,0.1^2,x0',0.05^2,2,0.9];
        LB = [eps,eps,0.00001,x0_LB,eps,eps,eps];
        UB = [1,0.3,0.3,x0_UB,0.3,20,1];
        if d == 4 %|| d == 9
            Xtrue = [0.9,0.05^2,0.1^2,x0',0.05^2,2,0.9];
        elseif d == 7 || d == 6
            Xtrue = [0.9,0.05^2,0.1^2,x0',0.05^2,2,0.99];
        end

        % V21
        Xtrue = [0.99,0.05^2,0.1^2,x0',0.05^2,0.01,0.9];
        LB = [eps,eps,0.00001,x0_LB,eps,eps,eps];
        UB = [1,0.3,0.3,x0_UB,0.3,0.5,1];

        % V22 % reducing lambda ranges
        Xtrue = [0.99,0.05^2,0.1^2,x0',0.05^2,0.01,0.9];
        LB = [eps,eps,0.00001,x0_LB,eps,eps,eps];
        UB = [1-eps,0.1,0.3,x0_UB,0.01,0.1,1-eps];

        % V23 % reducing lambda ranges
        Xtrue = [0.99,0.05^2,0.1^2,x0',0.05^2,0.01,0.999];
        LB = [0.85,eps,eps,x0_LB,eps,eps,0.9];
        UB = [1-eps,0.1,0.3,x0_UB,0.1,0.1,1-eps];

        % V24 % reducing lambda ranges
        Xtrue = [0.99,  0.01,   0.07,   x0',    2*eps,      0.01,   0.99];
        LB =    [0.85,  eps,    eps,    x0_LB,  eps,        eps,    0.9];
        UB =    [1-eps, 0.05,   0.1,    x0_UB,  0.0000001,  1,      1-eps];

        % V25 Like V20
        Xtrue = [0.99,0.05^2,0.1^2,x0',0.05^2,1,0.9];
        LB = [eps,eps,0.00001,x0_LB,eps,eps,eps];
        UB = [1,0.3,0.3,x0_UB,0.3,1,1];

        % V26 Like V20 --> seems ok not sim
        Xtrue = [0.99,0.05^2,0.1^2,x0',0.05^2,1,0.9];
        LB = [eps,eps,0.00001,x0_LB,eps,eps,eps];
        UB = [1,0.3,0.3,x0_UB,0.3,2,1];
        if d == 8 || d == 10
            Xtrue = [0.99,0.05^2,0.1^2,x0',0.05^2,1.5,0.9];
        end


        % V26 Like V20 --> seems ok not sim --> reducing pHAT increasing
        % sigmay range
        Xtrue = [0.99,0.05,0.1,x0',1e-8,1,0.9];
        LB = [eps,eps,0.00001,x0_LB,eps,eps,eps];
        UB = [1,0.3,0.5,x0_UB,1e-6,2,1];
        if d == 6
            Xtrue = [0.99,0.05,0.1,x0',1e-8,1.5,0.9];
        end

        % V27 Like V26 --> goood
        Xtrue = [0.99,0.05,0.1,x0',1e-8,.1,0.99];
        LB = [eps,eps,eps,x0_LB,eps,eps,eps];
        UB = [1,1e-5,5e-3,x0_UB,1e-5,1,1];
        if d == 6
            Xtrue = [0.99,0.05,0.1,x0',1e-8,1.5,0.9];
        end

        % V28 nope
        Xtrue = [0.999,5e-5,5e-3,x0',1e-6,.1,0.999];
        LB = [eps,eps,eps,x0_LB,eps,eps,eps];
        UB = [1,1e-3,1e-1,x0_UB,1e-4,1,1];
     
        % V29
        Xtrue = [0.999,5e-5,5e-3,x0',1e-6,.1,0.999];
        LB = [0.99,eps,eps,x0_LB,eps,eps,0.99];
        UB = [1,1e-5,1e-2,x0_UB,1e-5,1,1];

        % V30 very good also R^2
        Xtrue = [0.999,5e-6,8e-3,x0',1e-6,.01,0.999];
        LB = [0.99,eps,eps,x0_LB,eps,eps,0.99];
        UB = [1,1e-5,2e-2,x0_UB,1e-5,2,1];

         % V31 very good also R^2
        Xtrue = [0.999,5e-6,8e-3,x0',1e-6,.01,0.999];
        LB = [0.99,eps,eps,x0_LB,eps,eps,0.99];
        UB = [1,1e-5,4e-2,x0_UB,1e-5,3,1];

        % V32 very good also R^2
        Xtrue = [0.999,5e-6,8e-3,x0',1e-6,.01,0.999];
        LB = [0.99,eps,eps,x0_LB,eps,eps,0.99];
        UB = [1,1e-5,8e-2,x0_UB,1e-5,3,1];

        % V33 good also R^2
        Xtrue = [0.999,5e-6,3e-3,x0',1e-6,.01,0.999];
        LB = [0.99,eps,eps,x0_LB,eps,eps,0.99];
        UB = [1,1e-5,8e-2,x0_UB,1e-5,1,1];
        X0 = Xtrue;

        % V34  
        Xtrue = [0.999,5e-6,3e-3,x0',1e-6,.01,0.999];
        LB = [0.99,eps,eps,x0_LB,eps,eps,0.99];
        UB = [1,1e-5,8e-2,x0_UB,1e-5,0.1,1];
        X0 = Xtrue;

        % V35  
        Xtrue = [0.99,5e-6,3e-3,x0',1e-6,.01,0.99];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,1e-5,1e-1,x0_UB,1e-5,1,1-eps];
        

        % V36 nope
        Xtrue = [0.9,1e-6,1e-2,x0',1e-6,.1,0.9];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,2e-5,2e-2,x0_UB,2e-5,2,1-eps];
       

        % V37 good in sim global
        Xtrue = [0.9,1e-6,1e-2,x0',1e-6,.1,0.9];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,2e-5,5e-2,x0_UB,2e-5,2,1-eps];
        

        % V38 md nope
        Xtrue = [0.9,1e-6,1e-2,x0',1e-7,.1,0.9];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,2e-5,5e-2,x0_UB,1e-6,2,1-eps];

        % V39 no
        Xtrue = [0.99,1e-6,1e-2,x0',1e-7,.1,0.99];
        LB = [0.9,eps,eps,x0_LB,eps,eps,0.9];
        UB = [1-eps,2e-6,5e-2,x0_UB,1e-6,2,1-eps];

        % V40 no
        Xtrue = [0.99,1e-6,1e-2,x0',1e-7,.1,0.99];
        LB = [0.9,eps,eps,x0_LB,eps,eps,0.9];
        UB = [1-eps,2e-6,5e-2,x0_UB,1e-6,3,1-eps];

        % V41 no
        Xtrue = [0.99,1e-6,1e-2,x0',1e-7,.1,0.99];
        LB = [0.9,eps,eps,x0_LB,eps,eps,0.9];
        UB = [1-eps,2e-6,5e-2,x0_UB,1e-6,0.5,1-eps];

        % V42 from v37 
        Xtrue = [0.9,2e-6,5e-3,x0',2e-6,.1,0.9];
        LB = [0.9,1e-7,1e-7,x0_LB,1e-7,eps,0.9];
        UB = [0.9999,2e-5,5e-2,x0_UB,2e-5,2,0.9999];

        % V43 come 37
        Xtrue = [0.99,5e-6,3e-3,x0',1e-6,.01,0.99];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,2e-5,5e-2,x0_UB,2e-5,2,1-eps];

        % V44
        Xtrue = [0.99,5e-6,3e-3,x0',1e-6,.01,0.99];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [0.9999,2e-5,1e-1,x0_UB,2e-5,1,0.9999];

        % V45
        Xtrue = [0.99,5e-6,3e-3,x0',1e-6,.01,0.99];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,1e-5,1e-2,x0_UB,1e-5,0.5,1-eps];

        % V46 from simulations ---> not good 
        Xtrue = [0.99,5e-6,3e-3,x0',1e-6,.01,0.99];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [0.999,1e-6,1e-3,x0_UB,1e-7,2,0.9999];

        % V47 from simulations + x0 own action no
        % H: the other is doing my task
        for n = 1:nodes_num
            xmin_jeark(n) =  start(1) + (final(1)-start(1))*minimum_jerk_traj(tt(n+1)/duration);
            x0(n*2-1) =  xmin_jeark(n) ;
        end
        x01 = x0; x02 = x0;
        for n = 1:(nodes_num/2)
            ymin_jeark1(n) = start(2) + (VP1(2)-start(2))*minimum_jerk_traj(tt_half(n+1)/(duration/2));
            x01([n*2 (2*nodes_num-2*(n-1))]) = ymin_jeark1(n);
            ymin_jeark2(n) = start(2) + (VP2(2)-start(2))*minimum_jerk_traj(tt_half(n+1)/(duration/2));
            x02([n*2 (2*nodes_num-2*(n-1))]) = ymin_jeark2(n);
        end
        x01_LB = x01'-1e-4;
        x01_UB = x01'+1e-4;

        x02_LB = x02'-1e-4;
        x02_UB = x02'+1e-4;

        Xtrue1 = [0.99,1e-6,1e-2,x01',1e-6,.1,0.99];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,2e-5,5e-2,x01_UB,2e-5,2,1-eps];
        
        Xtrue2 =  [0.999,1e-6,1e-2,x02',1e-6,.1,0.99];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,2e-5,5e-2,x02_UB,2e-5,2,1-eps];
        
        % v48 good md no li
        Xtrue1 = [0.999,1e-6,5e-2,x01',1e-6,.1,0.999];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,2e-6,5e-2,x01_UB,2e-5,2,1-eps];
        
        Xtrue2 =  [0.999,1e-6,5e-2,x02',1e-6,.1,0.99];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,2e-6,5e-2,x02_UB,2e-5,2,1-eps];

        % v49 good li AND md --> E' LEI!
        Xtrue1 = [0.999,1e-6,5e-2,x01',1e-6,.1,0.999];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,2e-6,5e-2,x01_UB,2e-5,1.5,1-eps];
        
        Xtrue2 =  [0.999,1e-6,5e-2,x02',1e-6,.1,0.99];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,2e-6,5e-2,x02_UB,2e-5,1.5,1-eps];
       
        % v50 nopee
        Xtrue1 = [0.999,2e-6,5e-2,x01',1e-6,.1,0.999];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,2e-4,5e-2,x01_UB,2e-5,3,1-eps];
        
        Xtrue2 =  [0.999,2e-6,5e-2,x02',1e-6,.1,0.99];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,2e-4,5e-2,x02_UB,2e-5,3,1-eps];
       
        % v51 gni
        Xtrue1 = [0.999,1e-6,5e-2,x01',1e-6,.1,0.999];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,2e-6,5e-2,x01_UB,2e-5,1,1-eps];
        
        Xtrue2 =  [0.999,1e-6,5e-2,x02',1e-6,.1,0.99];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,2e-6,5e-2,x02_UB,2e-5,1,1-eps];

        % v52 no
        Xtrue1 = [0.999,1e-6,5e-2,x01',1e-6,.1,0.999];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,2e-5,5e-2,x01_UB,2e-5,1.5,1-eps];
        
        Xtrue2 =  [0.999,1e-6,5e-2,x02',1e-6,.1,0.99];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,2e-5,5e-2,x02_UB,2e-5,1.5,1-eps];

        opts = optimoptions('fmincon');
        opts.Display = 'iter';

        X01 = Xtrue1;
        X02 = Xtrue2;
        % partner 1
        X1 = fmincon(@(x) costfunction(x,task1,U2',U1'),X01,[],[],[],[],LB1,UB1,[],opts);
        % partner 2
        X2 = fmincon(@(x) costfunction(x,task2,U1',U2'),X02,[],[],[],[],LB1,UB1,[],opts);

        parnames_saving  = {'A','SigmaX','SigmaY','x0','P0','lambda0','a'};

        Pl1_Params = struct_params(parnames_saving,nodes_num,X1);
        Pl2_Params = struct_params(parnames_saving,nodes_num,X2);

        % Get cost
        parnames = {'A','B','C','D','SigmaX','SigmaY','x0','P0','lambda0','a'};
        pars1 = {X1(1)*eye(xsz1),[],eye(xsz1),eye(xsz1),X1(2)*eye(xsz1),X1(3)*eye(xsz1),X1(3+(1:xsz1)),X1(3+xsz1+1)*eye(xsz1),X1(3+xsz1+2),X1(3+xsz1+3)};
        model1 = playermodel(task1,pars1,parnames);
        [cost1,U1hat,Phat_U1,pi1] = pseudolik(model1,U2',U1');

        pars2 = {X2(1)*eye(xsz2),[],eye(xsz2),eye(xsz2),X2(2)*eye(xsz2),X2(3)*eye(xsz2),X2(3+(1:xsz2)),X2(3+xsz1+1)*eye(xsz2),X2(3+xsz2+2),X2(3+xsz2+3)};
        model2 = playermodel(task2,pars2,parnames);
        [cost2,U2hat,Phat_U2,pi2] = pseudolik(model2,U1',U2');

        % Visualization
        % Plotting weighted action
        for t = 1:size(pi1,1)
            uhat1_t(:,t) = zeros(nodes_num*2,1);
            uhat2_t(:,t) = zeros(nodes_num*2,1);

            for s = 1:size(pi1,2)
                uhat1_t(:,t) = uhat1_t(:,t) + pi1(t,s).*U1hat{t,s};
                uhat2_t(:,t) = uhat2_t(:,t) + pi2(t,s).*U2hat{t,s};
            end
        end

        figure
        subplot(121)
        plot_action(twoVPg,U1',U2')
        axis square
        title('u')
        subplot(122)
        plot_action(twoVPg,uhat1_t,uhat2_t)
        axis square
        title('uhat')
        set(gcf,'pos',[100 100 500 250])
        savefig(gcf,fullfile(datadir,['dyad' num2str(d) '_traj_' suffix '.fig']))

        % Plot Nodes and Y Nodes
        plot_nodes(twoVPg,nodes_num,U1,uhat1_t',U2,uhat2_t')
        savefig(gcf,fullfile(datadir,['Dyad' num2str(d) '_XYnodes_' suffix '.fig']))

        % Get Goodness of Simulations
        [dist1,dist2] = get_traj_differences(twoVPg,nodes_num,U1',uhat1_t,U2',uhat2_t);
        plot_traj_differences(twoVPg,nodes_num,U1',uhat1_t,U2',uhat2_t);

        [corrpl1,corrpl2,corrpl1x,corrpl1y,corrpl2x,corrpl2y] = get_corr2(nodes_num,U1',uhat1_t,U2',uhat2_t);

        % Saving
        save(fullfile(datadir, ['Dyads_' num2str(d) '_' suffix '.mat']),...
            "Pl1_Params","Pl2_Params","parnames_saving","tcE","tcL","cost1","cost2","U1hat","Phat_U1","pi1","U2hat","Phat_U2","pi2",...
            "dist1","dist2","corrpl1","corrpl2","corrpl1x","corrpl1y","corrpl2x","corrpl2y");
    end
    pause
end
%% Checking Parameters After fitting

if check
    GROUPcol = [0.7 0.7 0.3; 0.7 0.3 0.7; 0.3 0.7 0.7 ];
    sz = 30;
    Table_Params = {};
    Dyads_Num = [];
    % Loading Parameters
    for d = 1:n_dyads
        data_param = load(fullfile(datadir,['Dyads_' num2str(d) '_' suffix '.mat']));
        pl1_pos = (d-1)*2+1; pl2_pos = d*2;
        for p = 1:length(data_param.parnames_saving)
            param = data_param.parnames_saving{p}
            eval([param '1(d,:) = data_param.Pl1_Params.' param ';']);
            eval([param '2(d,:) = data_param.Pl2_Params.' param ';']);


            eval(['Table_Params{pl1_pos,p} = ' param '1(d,:);'])
            eval(['Table_Params{pl2_pos,p} = ' param '2(d,:);'])

        end

        Dyads_Num = [Dyads_Num; d; d];
    end
    Table_Params = cell2table(Table_Params);
    Table_Params.Properties.VariableNames = ["A","$\Sigma_x$","$\Sigma_y$","x0","$\hat{P_0}$","$\lambda$","$a_\lambda$"];
    Table_Params = removevars(Table_Params,["x0"]);    
    Group = {'H','H','H','H','H',...
        'H','H','H','H','H',...
        'VH','VH','VH','VH','VH',...
        'VH','VH','VH','VH','VH',...
        'PV','PV','PV','PV','PV',...
        'PV','PV','PV','PV','PV'}';

    Table_Params = addvars(Table_Params,Group,'Before',"A");
    Table_Params = addvars(Table_Params,Dyads_Num,'Before',"A");
    writetable(Table_Params,'ParametersAsymm_Table.xls')

    % BoxPlot for each parameter
    gr_labels = {'H','H','H','H','H',...
        'VH','VH','VH','VH','VH',...
        'PV','PV','PV','PV','PV',...
        'H','H','H','H','H',...
        'VH','VH','VH','VH','VH',...
        'PV','PV','PV','PV','PV'
        };
    params = {'A','SigmaX','SigmaY','lambda0','a','P0'};
    tit_params = {'A','$\Sigma_x$','$\Sigma_y$','$\lambda$','a','Phat'};

    params = {'SigmaY'};
    tit_params = {'$\Sigma_y$'};

    for pp = 1:length(params)

        eval(['p1 =' params{pp} num2str(1) ';']);
        eval(['p2 =' params{pp} num2str(2) ';']);

        eval(['y = [' params{pp} '1; ' params{pp} '2];'])

        % stat
        [p,tbl,stats] = anova1(y,gr_labels);
        title(params{pp})
        xticklabels({'H','VH','PV'})
        [results,m,h] = multcompare(stats);
        [results_corr,m_corr,h_corr] = multcompare(stats,'CriticalValueType','hsd');
        tbl = array2table(results,"VariableNames", ...
            ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

        % boxplot
        maxval = max([p1; p2]);
        minval = min([p1; p2]);

        figure
        set(gcf,'pos',[0 0 160 280])
        boxplot([reshape(p1,n_dyads_group,n_groups);reshape(p2,n_dyads_group,n_groups)])

        h = findobj(gca,'Tag','Box');
        mm = findobj(gca,'Tag','Median');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),GROUPcol(j,:),'FaceAlpha',.3,...
                'Edgecolor',GROUPcol(j,:));
        end

        hold on
        %         plot(1:n_groups,reshape(p1,n_dyads,n_groups)','k*')
        %         plot(1:n_groups,reshape(p2,n_dyads,n_groups)','k*')
        hold on
        plot(ones(1,10)+(rand(size(ones(1,10)))-0.5).*0.5,[p1(1:5); p2(1:5)],'.','color',GROUPcol(3,:),'markersize',10)
        plot(2.*ones(1,10)+(rand(size(ones(1,10)))-0.5).*0.5,[p1(6:10); p2(6:10)],'.','color',GROUPcol(2,:),'markersize',10)
        plot(3.*ones(1,10)+(rand(size(ones(1,10)))-0.5).*0.5,[p1(11:15); p2(11:15)],'.','color',GROUPcol(1,:),'markersize',10)

        box off
        xticklabels({'H','VH','PV'})
        %title(['p = ' num2str(p)])
        ylabel(tit_params{pp},'interpreter','latex')
        ylim([minval-0.01*minval maxval+0.01*minval])
        set(gca, 'FontName', 'times');
    end

    % Figure SigmaY-Lambda
    % figure
    % hold on
    % scatter([SigmaY1(subH) SigmaY2(subH)],[lambda01(subH) lambda02(subH)],sz,GROUPcol(3,:),'filled')
    % scatter([SigmaY1(subVH) SigmaY2(subVH)],[lambda01(subVH) lambda02(subVH)],sz,GROUPcol(2,:),'filled')
    % scatter([SigmaY1(subPV) SigmaY2(subPV)],[lambda01(subPV) lambda02(subPV)],sz,GROUPcol(1,:),'filled')
    % ylabel('Lambda')
    % xlabel('Sigma_y')
    % set(gca, 'FontName', 'times');
    % set(gcf,'pos',[0 0 300 300])

    save("EstSigmaY.mat",'SigmaY1','SigmaY2')
end

pause
%% Simulation
if simulate_single_dyads

    clear corrpl1 corrpl2
    load(fullfile(datadir, ['forces.mat']));
    indicators = {'md11','md22','md12','md21','tc11','tc22','te1','te2','ts1','ts2'};
    unconnected = find(forces(end,:) == 0);

    for i = 1:length(indicators)
        load(fullfile(datadir, [indicators{i} '.mat']));
        eval([indicators{i} '= ind_ts;' ]);
        eval([indicators{i} '(:,unconnected)= [];' ]);
    end

    for d = 1:n_dyads
        % Defining task specs
        [tcE,tcL, duration, yrange, start, final, VP1, VP2, tgterror,VPerror] = get_2VPtask_specs('Asymm');

        % Reading action data
        if nnodes_dataset== 5
            data_action = readtable(fullfile(datadir,['dyad' num2str(d) '_5nodes.dat']));
        else
            data_action = readtable(fullfile(datadir,['dyad' num2str(d) '.dat']));
        end

        % Removing unconnected trials
        %unconnected = find(forces(d,:) == 0);
        data_action(unconnected,:) = [];

        nodes_num = size(data_action,2)/4; % number of nodes in spline approximation
        T = size(data_action,1); % number of plays

        % Partner 1 actions
        U1 = table2array(data_action(:,1:(2*nodes_num)));
        % Partner 2 actions
        U2 = table2array(data_action(:,2*nodes_num+(1:(2*nodes_num))));

        if nnodes_dataset == 5
            U1(:,1:2:9) = U1(:,1:2:9)-0.05;
            U2(:,1:2:9) = U2(:,1:2:9)-0.05;
        end

        twoVPg = twoVP_game(nodes_num,VP1,VP2,tcE,tcL,start,final,duration,VPerror,tgterror);
        task1 = get(twoVPg,'task1');
        task2 = get(twoVPg,'task2');

        xsz1 = get(task1,'xsize');
        xsz2 = get(task2,'xsize');

        % DISPLAY Nash equilibria
        game = quadraticgame(task1,task2);
        [u1ne,u2ne] = nash_equilibrium(game);
        ne_label = {'VP_1 \rightarrow VP_2','VP_2 \rightarrow VP_1','NC'};

        for n = 1:length(u1ne)
            n
             figure('pos',[100 100 400 400])
             plot_action(twoVPg,u1ne{n,n},u2ne{n,n})

            J1 = cost(task1,u1ne{n,n},u2ne{n,n});
            J2 = cost(task2,u2ne{n,n},u1ne{n,n});

             J1{n}
             J2{n}
             title(sprintf('Nash equilibrium %d: %s',n,ne_label{n}));
             % text(0,0.05,num2str(J1{n}),'color','b');
             % text(0,0.04,num2str(J2{n}),'color','r');
             ylim([-0.06 0.06])
             xlim([-0.06 0.06])
             set(gcf,'Position',[100 100 300 300])
        end

        % Initialize action selection variables
        collab{d} = zeros(n_sim,T);collab12{d} = zeros(n_sim,T);collab21{d} = zeros(n_sim,T);
        oscil{d} = zeros(n_sim,T);oscil12_21{d} = zeros(n_sim,T);oscil21_12{d} = zeros(n_sim,T);
        bothignore{d} = zeros(n_sim,T);oneignore{d} = zeros(n_sim,T); nocollab{d} = zeros(n_sim,T);

        % Define the model
        % Loading Parameters
        xsz = nodes_num*2;
        data_param = load(fullfile(datadir,['Dyads_' num2str(d) '_' suffix '.mat']));
        [H1,L1,sigmay1,A1,sigmax1,P01,x01,lambda1,decay_rate1] = define_player_model_from_fitting(data_param.Pl1_Params,xsz,duration,nodes_num,final,start,VP1);
        [H2,L2,sigmay2,A2,sigmax2,P02,x02,lambda2,decay_rate2] = define_player_model_from_fitting(data_param.Pl2_Params,xsz,duration,nodes_num,final,start,VP2);


        % Observer INITIALIZATION
        % Initialize the observer(s)
        observer1 = observer(A1,[],H1,L1,sigmax1,sigmay1,x01,P01);
        observer2 = observer(A2,[],H2,L2,sigmax2,sigmay2,x02,P02);

        % Initialize the simulation
        P1_prior = cell(T,1);
        P2_prior = cell(T,1);
        x1_prior = zeros(xsz1,T);
        x2_prior = zeros(xsz2,T);

        x1_prior(:,1) = x01;
        P1_prior{1} = P01;

        x2_prior(:,1) = x02;
        P2_prior{1} = P02;

        lambda1(1) = lambda1;
        lambda2(1) = lambda2;


        for s = 1:n_sim

            contr1 = controller(task1,lambda1(1),decay_rate1);
            contr2 = controller(task2,lambda2(1),decay_rate2);

            for t=1:T
                t
                % Generate new actions
                switch gen_action_mode
                    case 'stochastic'
                        [u1(:,t),decision1(:,t), contr1] = generate_action(contr1, x1_prior(:,t),P1_prior{t},lambda1(t));
                        [u2(:,t),decision2(:,t), contr2] = generate_action(contr2, x2_prior(:,t),P2_prior{t},lambda2(t));
                    case 'deterministic'
                        [u1(:,t),decision1(:,t), contr1] = generate_action_deterministic(contr1, x1_prior(:,t),P1_prior{t},lambda1(t));
                        [u2(:,t),decision2(:,t), contr2] = generate_action_deterministic(contr2, x2_prior(:,t),P2_prior{t},lambda2(t));
                end
                prior1(t,:) = get(contr1,'prior');
                prior2(t,:) = get(contr2,'prior');

                % Generate sensory feedback
                switch mode
                    case 'single'
                        y1(:,t) = generate_sensory(observer1,U2(t,:)',u1(:,t),'sim');
                        y2(:,t) = generate_sensory(observer2,U1(t,:)',u2(:,t),'sim');
                    case 'dyad'
                        y1(:,t) = generate_sensory(observer1,u2(:,t),u1(:,t),'sim');
                        y2(:,t) = generate_sensory(observer2,u1(:,t),u2(:,t),'sim');
                end
                % Update state observers: dynamic
                [x1_post(:,t),P1_post{t},x1_prior(:,t+1),P1_prior{t+1},Kalman1{t}] = kalman_estimate(observer1,x1_prior(:,t),P1_prior{t},y1(:,t),u1(:,t));
                [x2_post(:,t),P2_post{t},x2_prior(:,t+1),P2_prior{t+1},Kalman2{t}] = kalman_estimate(observer2,x2_prior(:,t),P2_prior{t},y2(:,t),u2(:,t));

                % Update temperature
                lambda1(t+1) = decay_rate1*lambda1(t);
                lambda2(t+1) = decay_rate2*lambda2(t);
            end
            u1sim{d}(:,:,s) = u1;
            u2sim{d}(:,:,s) = u2;

            %save('u1','u2','prior1','prior2','x1_prior','x2_prior','y1','y2',fullfile(savedir,fname));

            % SAVE simulation data
            data = table(u1',u2','VariableNames',{'u1','u2'});
            %             writetable(data,fullfile(savedir,fname));

            % plot actions
            figure
            plot_action(twoVPg,u1(:,(end-9):end),u2(:,(end-9):end));
            set(gcf,'pos',[0 500 175 175])

            % figure
            % subplot(121)
            % plot_action(twoVPg,U1',U2');
            % subplot(122)
            % plot_action(twoVPg,u1,u2);
            % set(gcf,'pos',[0 500 300 150])


           %plot_animation(twoVPg,u1,u2);

            % plot speed profile
            %plot_speed(twoVPg,u1,u2);

            % plot acceleration profile
            %plot_acceleration(twoVPg,u1,u2);

            % plot minimum distance
            %plot_minimum_distance(twoVPg,u1,u2);
            [md11_sim{d}(s,:),md12_sim{d}(s,:),md21_sim{d}(s,:),md22_sim{d}(s,:),...
                tc11_sim{d}(s,:),tc12_sim{d}(s,:),tc21_sim{d}(s,:),tc22_sim{d}(s,:)] = get_minimum_distance(twoVPg,u1,u2);

            if s == 1
                exp_sig = [md11(d,:);md12(d,:);md21(d,:);md22(d,:)];
                sim_sig = [md11_sim{d}(s,:);md12_sim{d}(s,:);md21_sim{d}(s,:);md22_sim{d}(s,:)];
                %plot_md_exp_sim(twoVPg,exp_sig,sim_sig,1:T,[0 0.05],{'MD_{11}','MD_{12}','MD_{21}','MD{22}'})
            end
            %saveas(gcf,fullfile(savepath,'MD.png'))

            % Get Leadership Index
            [LI_11{d}(s,:), LI_12{d}(s,:), LI_21{d}(s,:), LI_22{d}(s,:)] = get_LI(twoVPg,u1,u2,tcE*duration,tcL*duration,duration);
            %plot_LI(twoVPg,u1,u2,tcE*duration,tcL*duration,duration);

            % Get Spatial Variability
            block_trials=10;
            [spvar1{d}(s,:),spvar2{d}(s,:)] = get_spvar(twoVPg,u1,u2,block_trials);

            % Get Goodness of Simulations
            [dist1_sim{d}(s,:),dist2_sim{d}(s,:)] = get_traj_differences(twoVPg,nodes_num,U1',u1,U2',u2);
            %plot_traj_differences(twoVPg,nodes_num,U1',u1,U2',u2);

            [corrpl1{d}(s,:),corrpl2{d}(s,:),corrpl1x{d}(s,:),corrpl1y{d}(s,:),corrpl2x{d}(s,:),corrpl2y{d}(s,:)] = get_corr2(nodes_num,U1',u1,U2',u2);

            u1temp = u1sim{d};
            u2temp = u2sim{d};
            dist1temp = dist1_sim{d}; dist2temp = dist2_sim{d};
            spvar1temp = spvar1{d}; spvar2temp = spvar2{d};
            corrcoeff1temp = corrpl1{d}; corrcoeff1xtemp = corrpl1x{d}; corrcoeff1ytemp = corrpl1y{d}; 
            corrcoeff2temp = corrpl2{d}; corrcoeff2xtemp = corrpl2x{d}; corrcoeff2ytemp = corrpl2y{d}; 
            mindist12 = md12_sim{d}; mindist21 = md21_sim{d};
            save(fullfile(rootfolder,'data','simulations','VPAsymm',['dyad_' num2str(d) '_V' num2str(v) '_' mode '.mat']),...
                'U1','U2','u1temp','u2temp',"dist1temp","dist2temp", "spvar2temp","spvar2temp","mindist12","mindist21",...
                "corrcoeff1temp",'corrcoeff2temp',"corrcoeff1xtemp",'corrcoeff2xtemp',"corrcoeff1ytemp",'corrcoeff2ytemp')

        end
        %close all
    end
    

    %% Plot Group MDs
    colorgroup = [0.3 0.7 0.7; 0.7 0.3 0.7; 0.7 0.7 0.3 ];
    alpha = .3
    indicators = {'md12_sim','md21_sim'}
    gg = {subH,subVH,subPV};
    figure
    for i = 1:length(indicators)
        eval(['curr_ind =' indicators{i} ';'])

        subplot(1,2,i)
        hold on
        n_sub = n_dyads_group;
        for g = 1:n_groups
            sub_Group = gg{g};
            ind_group{g} = [];
            for ii = sub_Group
                ind_group{g} = [ind_group{g}; curr_ind{ii}];
            end
        end
        for g = 1:n_groups
            patch([1:T T:-1:1],[mean(ind_group{g},1)-std(ind_group{g},1)./sqrt(n_sub) fliplr(mean(ind_group{g},1)+std(ind_group{g},1)./sqrt(n_sub))],...
                colorgroup(g,:),'facealpha',alpha,'edgecolor',colorgroup(g,:),'edgealpha',alpha)
        end
        for g = 1:n_groups
            s(g) = plot(mean(ind_group{g}),'linewidth',1.5,'color',colorgroup(g,:))
        end


        ylim([0 0.05])
        xlim([0 T+1])
        xlabel('Trial')
        ylabel('distance [m]')
        if i == length(indicators)
            subplot(1,2,i)
            legend(s,{'H','VH','PV'})
            legend boxoff
            set(gca,'FontName','Times')

        end

    end
    set(gcf,'pos',[0 0 360 250])

    %% Plot Group MDs, LI
    colorgroup = [0.3 0.7 0.7; 0.7 0.3 0.7; 0.7 0.7 0.3 ];
    alpha = .3
    indicators = {'md12_sim','md21_sim','LI_11','LI_12','LI_21','LI_22'};
    yl = [0 0.05;0 0.05; -0.5 0.5;-0.5 0.5;-0.5 0.5;-0.5 0.5];
    ylab = {'distance [m]','distance [m]','[w]','[w]','[w]','[w]'};
    tit = {'MD12','MD21','LI_{11}','LI_{12}','LI_{21}','LI_{22}'};
    gg = {subH,subVH,subPV};
    
    for i = 1:length(indicators)
        eval(['curr_ind =' indicators{i} ';'])
          
        for d = 1:n_dyads
            curr_ind{d} = movmean(curr_ind{d},1);
        end
 
        figure
        %subplot(1,2,i)
        hold on
        n_sub = n_dyads_group;
        for g = 1:n_groups
            sub_Group = gg{g};
            ind_group{g} = [];
            for ii = sub_Group
                ind_group{g} = [ind_group{g}; curr_ind{ii}];
            end
        end
        for g = 1:n_groups
            patch([1:T T:-1:1],[mean(ind_group{g},1)-std(ind_group{g},1)./sqrt(n_sub) fliplr(mean(ind_group{g},1)+std(ind_group{g},1)./sqrt(n_sub))],...
                colorgroup(g,:),'facealpha',alpha,'edgecolor',colorgroup(g,:),'edgealpha',alpha)
        end
        for g = 1:n_groups
            s(g) = plot(mean(ind_group{g}),'linewidth',1.5,'color',colorgroup(g,:))
        end

        ylim(yl(i,:))
        xlim([0 T+1])
        xlabel('Trial')
        ylabel(ylab{i})
        title(tit{i})
        %if i == length(indicators)
            %subplot(1,2,i)
            legend(s,{'H','VH','PV'})
            legend boxoff
        %end
        set(gca,'fontname','Times')
        set(gcf,'pos',[0 0 160 210])
        set(gcf,'pos',[0 0 125 204])
        set(gcf,'pos',[0 0 170 250])
    set(gcf,'pos',[0 0 150 210])

    end

    % Leadership bar
    for d = 1:n_dyads
        LI11(d,:) = LI_11{d};
        LI12(d,:) = LI_12{d};
        LI21(d,:) = LI_21{d};
        LI22(d,:) = LI_22{d};
        LI_VP1(d,:) = LI21(d,:)- LI11(d,:);
        LI_VP2(d,:) = LI12(d,:) - LI22(d,:);
        LI_pl1(d,:) = LI12(d,:)- LI11(d,:);
        LI_pl2(d,:) = LI21(d,:) - LI22(d,:);
    end

    for d = 1:n_dyads
        LI_11_av(d) = mean( LI11(d,(end-9):end));
        LI_12_av(d) = mean( LI12(d,(end-9):end));
        LI_21_av(d) = mean( LI21(d,(end-9):end));
        LI_22_av(d) = mean( LI22(d,(end-9):end));
        LI_VP1_av(d) = mean( LI_VP1(d,(end-9):end));
        LI_VP2_av(d) = mean( LI_VP2(d,(end-9):end));
        LI_pl1_av(d) = mean( LI_pl1(d,(end-9):end));
        LI_pl2_av(d) = mean( LI_pl2(d,(end-9):end));
    end


inds = {'LI_11_av','LI_12_av','LI_21_av','LI_22_av'};%,'LI_VP1_av','LI_VP2_av','LI_pl1_av','LI_pl2_av'};
titles = {'LI_{11}','LI_{12}','LI_{21}','LI_{22}','\Delta LI_1','\Delta LI_2','\Delta LI_1','\Delta LI_2'};
ylims = [-0.5 0.5; -0.5 0.5; -0.5 0.5; -0.5 0.5];% -0.2 0.6; -0.2 0.6; -0.-0.2 0.6; -0.2 0.6];

for ii = 1:length(inds)
    eval(['LI =' inds{ii} ';']);

    figure
    title(titles{ii})
    hold on
    for g = 1:n_groups
        curr_dyads = (g-1)*5+(1:5);
        errorbar(g,mean(LI(curr_dyads)),std(LI(curr_dyads))./sqrt(5),'k')
        bar(g,mean(LI(curr_dyads)),'Facecolor',colorgroup(g,:))
    end
    %sigline([1 3],[],[],[])

    ylabel('[w]')
    xticks(1:3)
    xticklabels({'H','VH','PV'})
    ylim(ylims(ii,:))
    set(gca,'fontname','times')
    set(gcf,'pos',[0 0 150 210])

    [p,tbl,stats] = anova1(reshape(LI,5,3))
    c = multcompare(stats)
end


end

%% Simulating average parameters

if simulate_global_groups 


    clear corrpl1 corrpl2
    load(fullfile(datadir, ['forces.mat']));
    indicators = {'md11','md22','md12','md21','tc11','tc22','te1','te2','ts1','ts2'};
    unconnected = find(forces(end,:) == 0);

    for i = 1:length(indicators)
        load(fullfile(datadir, [indicators{i} '.mat']));
        eval([indicators{i} '= ind_ts;' ]);
        eval([indicators{i} '(:,unconnected)= [];' ]);
    end

    Params = defining_player_model_from_fitting_bygroup(datadir,suffix,n_dyads,n_groups,n_dyads_group)

    for d = 1:n_dyads
        % Defining task specs
        [tcE,tcL, duration, yrange, start, final, VP1, VP2, tgterror,VPerror] = get_2VPtask_specs('Asymm');

        % Reading action data
        if nnodes_dataset== 5
            data_action = readtable(fullfile(datadir,['dyad' num2str(d) '_5nodes.dat']));
        else
            data_action = readtable(fullfile(datadir,['dyad' num2str(d) '.dat']));
        end

        % Removing unconnected trials
        %unconnected = find(forces(d,:) == 0);
        data_action(unconnected,:) = [];

        nodes_num = size(data_action,2)/4; % number of nodes in spline approximation
        T = size(data_action,1); % number of plays

        % Partner 1 actions
        U1 = table2array(data_action(:,1:(2*nodes_num)));
        % Partner 2 actions
        U2 = table2array(data_action(:,2*nodes_num+(1:(2*nodes_num))));

        if nnodes_dataset == 5
            U1(:,1:2:9) = U1(:,1:2:9)-0.05;
            U2(:,1:2:9) = U2(:,1:2:9)-0.05;
        end

        twoVPg = twoVP_game(nodes_num,VP1,VP2,tcE,tcL,start,final,duration,VPerror,tgterror);
        task1 = get(twoVPg,'task1');
        task2 = get(twoVPg,'task2');

        xsz1 = get(task1,'xsize');
        xsz2 = get(task2,'xsize');

        % DISPLAY Nash equilibria
        game = quadraticgame(task1,task2);
        [u1ne,u2ne] = nash_equilibrium(game);
        ne_label = {'VP_1 \rightarrow VP_2','VP_2 \rightarrow VP_1','NC'};

        for n = 1:length(u1ne)
            n
            % figure('pos',[100 100 400 400])
            % plot_action(twoVPg,u1ne{n,n},u2ne{n,n})

            J1 = cost(task1,u1ne{n,n},u2ne{n,n});
            J2 = cost(task2,u2ne{n,n},u1ne{n,n});

            % J1{n}
            % J2{n}
            % title(sprintf('Nash equilibrium %d: %s',n,ne_label{n}));
            % text(0,0.05,num2str(J1{n}),'color','b');
            % text(0,0.04,num2str(J2{n}),'color','r');
            % ylim([-0.06 0.06])
            % xlim([-0.06 0.06])
            % set(gcf,'Position',[100 100 300 300])
        end

        % Initialize action selection variables
        collab{d} = zeros(n_sim,T);collab12{d} = zeros(n_sim,T);collab21{d} = zeros(n_sim,T);
        oscil{d} = zeros(n_sim,T);oscil12_21{d} = zeros(n_sim,T);oscil21_12{d} = zeros(n_sim,T);
        bothignore{d} = zeros(n_sim,T);oneignore{d} = zeros(n_sim,T); nocollab{d} = zeros(n_sim,T);

        % Define the model
        % Loading Parameters
        xsz = nodes_num*2;
        %data_param = load(fullfile(datadir,['Dyads_' num2str(d) '_' suffix '.mat']));
        [H1,L1,sigmay1,A1,sigmax1,P01,x01,lambda1,decay_rate1] = define_player_model_from_fitting(Params{d},xsz,duration,nodes_num,final,start,VP1);
        [H2,L2,sigmay2,A2,sigmax2,P02,x02,lambda2,decay_rate2] = define_player_model_from_fitting(Params{d},xsz,duration,nodes_num,final,start,VP2);


        % Observer INITIALIZATION
        % Initialize the observer(s)
        observer1 = observer(A1,[],H1,L1,sigmax1,sigmay1,x01,P01);
        observer2 = observer(A2,[],H2,L2,sigmax2,sigmay2,x02,P02);

        % Initialize the simulation
        P1_prior = cell(T,1);
        P2_prior = cell(T,1);
        x1_prior = zeros(xsz1,T);
        x2_prior = zeros(xsz2,T);

        x1_prior(:,1) = x01;
        P1_prior{1} = P01;

        x2_prior(:,1) = x02;
        P2_prior{1} = P02;

        lambda1(1) = lambda1;
        lambda2(1) = lambda2;


        for s = 1:n_sim

            contr1 = controller(task1,lambda1(1),decay_rate1);
            contr2 = controller(task2,lambda2(1),decay_rate2);

            for t=1:T
                t
                % Generate new actions
                switch gen_action_mode
                    case 'stochastic'
                        [u1(:,t),decision1(:,t), contr1] = generate_action(contr1, x1_prior(:,t),P1_prior{t},lambda1(t));
                        [u2(:,t),decision2(:,t), contr2] = generate_action(contr2, x2_prior(:,t),P2_prior{t},lambda2(t));
                    case 'deterministic'
                        [u1(:,t),decision1(:,t), contr1] = generate_action_deterministic(contr1, x1_prior(:,t),P1_prior{t},lambda1(t));
                        [u2(:,t),decision2(:,t), contr2] = generate_action_deterministic(contr2, x2_prior(:,t),P2_prior{t},lambda2(t));
                end
                prior1(t,:) = get(contr1,'prior');
                prior2(t,:) = get(contr2,'prior');

                % Generate sensory feedback
                switch mode
                    case 'single'
                        y1(:,t) = generate_sensory(observer1,U2(t,:)',u1(:,t),'sim');
                        y2(:,t) = generate_sensory(observer2,U1(t,:)',u2(:,t),'sim');
                    case 'dyad'
                        y1(:,t) = generate_sensory(observer1,u2(:,t),u1(:,t),'sim');
                        y2(:,t) = generate_sensory(observer2,u1(:,t),u2(:,t),'sim');
                end
                % Update state observers: dynamic
                [x1_post(:,t),P1_post{t},x1_prior(:,t+1),P1_prior{t+1},Kalman1{t}] = kalman_estimate(observer1,x1_prior(:,t),P1_prior{t},y1(:,t),u1(:,t));
                [x2_post(:,t),P2_post{t},x2_prior(:,t+1),P2_prior{t+1},Kalman2{t}] = kalman_estimate(observer2,x2_prior(:,t),P2_prior{t},y2(:,t),u2(:,t));

                % Update temperature
                lambda1(t+1) = decay_rate1*lambda1(t);
                lambda2(t+1) = decay_rate2*lambda2(t);
            end
            u1sim{d}(:,:,s) = u1;
            u2sim{d}(:,:,s) = u2;

            %save('u1','u2','prior1','prior2','x1_prior','x2_prior','y1','y2',fullfile(savedir,fname));

            % SAVE simulation data
            data = table(u1',u2','VariableNames',{'u1','u2'});
            %             writetable(data,fullfile(savedir,fname));

            % plot actions
            figure
            plot_action(twoVPg,u1(:,(end-9):end),u2(:,(end-9):end));
            set(gcf,'pos',[0 500 150 150])
            % figure
            % subplot(121)
            % plot_action(twoVPg,U1',U2');
            % subplot(122)
            % plot_action(twoVPg,u1,u2);
            % set(gcf,'pos',[0 500 300 150])


           %plot_animation(twoVPg,u1,u2);

            % plot speed profile
            %plot_speed(twoVPg,u1,u2);

            % plot acceleration profile
            %plot_acceleration(twoVPg,u1,u2);

            % plot minimum distance
            %plot_minimum_distance(twoVPg,u1,u2);
            [md11_sim{d}(s,:),md12_sim{d}(s,:),md21_sim{d}(s,:),md22_sim{d}(s,:),...
                tc11_sim{d}(s,:),tc12_sim{d}(s,:),tc21_sim{d}(s,:),tc22_sim{d}(s,:)] = get_minimum_distance(twoVPg,u1,u2);

            if s == 1
                exp_sig = [md11(d,:);md12(d,:);md21(d,:);md22(d,:)];
                sim_sig = [md11_sim{d}(s,:);md12_sim{d}(s,:);md21_sim{d}(s,:);md22_sim{d}(s,:)];
                %plot_md_exp_sim(twoVPg,exp_sig,sim_sig,1:T,[0 0.05],{'MD_{11}','MD_{12}','MD_{21}','MD{22}'})
            end
            %saveas(gcf,fullfile(savepath,'MD.png'))

            % Get Leadership Index
            [LI_11{d}(s,:), LI_12{d}(s,:), LI_21{d}(s,:), LI_22{d}(s,:)] = get_LI(twoVPg,u1,u2,tcE*duration,tcL*duration,duration);
            %plot_LI(twoVPg,u1,u2,tcE*duration,tcL*duration,duration);

            % Get Spatial Variability
            block_trials=10;
            [spvar1{d}(s,:),spvar2{d}(s,:)] = get_spvar(twoVPg,u1,u2,block_trials);

            % Get Goodness of Simulations
            [dist1_sim{d}(s,:),dist2_sim{d}(s,:)] = get_traj_differences(twoVPg,nodes_num,U1',u1,U2',u2);
            %plot_traj_differences(twoVPg,nodes_num,U1',u1,U2',u2);

            [corrpl1{d}(s,:),corrpl2{d}(s,:),corrpl1x{d}(s,:),corrpl1y{d}(s,:),corrpl2x{d}(s,:),corrpl2y{d}(s,:)] = get_corr2(nodes_num,U1',u1,U2',u2);

            u1temp = u1sim{d};
            u2temp = u2sim{d};
            dist1temp = dist1_sim{d}; dist2temp = dist2_sim{d};
            spvar1temp = spvar1{d}; spvar2temp = spvar2{d};
            corrcoeff1temp = corrpl1{d}; corrcoeff1xtemp = corrpl1x{d}; corrcoeff1ytemp = corrpl1y{d}; 
            corrcoeff2temp = corrpl2{d}; corrcoeff2xtemp = corrpl2x{d}; corrcoeff2ytemp = corrpl2y{d}; 
            mindist12 = md12_sim{d}; mindist21 = md21_sim{d};
            save(fullfile(rootfolder,'data','simulations','VPAsymm',['dyad_' num2str(d) '_V' num2str(v) '_' mode '.mat']),...
                'U1','U2','u1temp','u2temp',"dist1temp","dist2temp", "spvar2temp","spvar2temp","mindist12","mindist21",...
                "corrcoeff1temp",'corrcoeff2temp',"corrcoeff1xtemp",'corrcoeff2xtemp',"corrcoeff1ytemp",'corrcoeff2ytemp')

        end
        %close all
    end
    
    %% Plot Group MDs, LI
    colorgroup = [0.3 0.7 0.7; 0.7 0.3 0.7; 0.7 0.7 0.3 ];
    alpha = .3
    indicators = {'md12_sim','md21_sim','LI_11','LI_12','LI_21','LI_22'};
    yl = [0 0.05;0 0.05; -0.5 0.5;-0.5 0.5;-0.5 0.5;-0.5 0.5];
    ylab = {'distance [m]','distance [m]','[w]','[w]','[w]','[w]'};
    tit = {'MD12','MD21','LI_{11}','LI_{12}','LI_{21}','LI_{22}'};
    gg = {subH,subVH,subPV};
    
    for i = 1:length(indicators)
        eval(['curr_ind =' indicators{i} ';'])
          
        for d = 1:n_dyads
            curr_ind{d} = movmean(curr_ind{d},1);
        end
 
        figure
        %subplot(1,2,i)
        hold on
        n_sub = n_dyads_group;
        for g = 1:n_groups
            sub_Group = gg{g};
            ind_group{g} = [];
            for ii = sub_Group
                ind_group{g} = [ind_group{g}; curr_ind{ii}];
            end
        end
        for g = 1:n_groups
            patch([1:T T:-1:1],[mean(ind_group{g},1)-std(ind_group{g},1)./sqrt(n_sub) fliplr(mean(ind_group{g},1)+std(ind_group{g},1)./sqrt(n_sub))],...
                colorgroup(g,:),'facealpha',alpha,'edgecolor',colorgroup(g,:),'edgealpha',alpha)
        end
        for g = 1:n_groups
            s(g) = plot(mean(ind_group{g}),'linewidth',1.5,'color',colorgroup(g,:))
        end

        ylim(yl(i,:))
        xlim([0 T+1])
        xlabel('Trial')
        ylabel(ylab{i})
        title(tit{i})
        %if i == length(indicators)
            %subplot(1,2,i)
            legend(s,{'H','VH','PV'})
            legend boxoff
        %end
        set(gca,'fontname','Times')
        set(gcf,'pos',[0 0 160 210])
        set(gcf,'pos',[0 0 125 204])
        set(gcf,'pos',[0 0 170 250])

    end

    % Leadership bar
    for d = 1:n_dyads
        LI11(d,:) = LI_11{d};
        LI12(d,:) = LI_12{d};
        LI21(d,:) = LI_21{d};
        LI22(d,:) = LI_22{d};
        LI_VP1(d,:) = LI21(d,:)- LI11(d,:);
        LI_VP2(d,:) = LI12(d,:) - LI22(d,:);
        LI_pl1(d,:) = LI12(d,:)- LI11(d,:);
        LI_pl2(d,:) = LI21(d,:) - LI22(d,:);
    end

    for d = 1:n_dyads
        LI_11_av(d) = mean( LI11(d,(end-9):end));
        LI_12_av(d) = mean( LI12(d,(end-9):end));
        LI_21_av(d) = mean( LI21(d,(end-9):end));
        LI_22_av(d) = mean( LI22(d,(end-9):end));
        LI_VP1_av(d) = mean( LI_VP1(d,(end-9):end));
        LI_VP2_av(d) = mean( LI_VP2(d,(end-9):end));
        LI_pl1_av(d) = mean( LI_pl1(d,(end-10):end));
        LI_pl2_av(d) = mean( LI_pl2(d,(end-10):end));
    end


inds = {'LI_11_av','LI_12_av','LI_21_av','LI_22_av'};%,'LI_VP1_av','LI_VP2_av','LI_pl1_av','LI_pl2_av'};
titles = {'LI_{11}','LI_{12}','LI_{21}','LI_{22}','\Delta LI_1','\Delta LI_2','\Delta LI_1','\Delta LI_2'};
ylims = [-0.5 0.5; -0.5 0.5; -0.5 0.5; -0.5 0.5];% -0.2 0.6; -0.2 0.6; -0.-0.2 0.6; -0.2 0.6];

for ii = 1:length(inds)
    eval(['LI =' inds{ii} ';']);

    figure
    title(titles{ii})
    hold on
    for g = 1:n_groups
        curr_dyads = (g-1)*5+(1:5);
        errorbar(g,mean(LI(curr_dyads)),std(LI(curr_dyads))./sqrt(5),'k')
        bar(g,mean(LI(curr_dyads)),'Facecolor',colorgroup(g,:))
    end
    %sigline([1 3],[],[],[])

    ylabel('[w]')
    xticks(1:3)
    xticklabels({'H','VH','PV'})
    ylim(ylims(ii,:))
    set(gca,'fontname','times')
    set(gcf,'pos',[0 0 125 204])

    [p,tbl,stats] = anova1(reshape(LI,5,3))
    c = multcompare(stats)
end

end





