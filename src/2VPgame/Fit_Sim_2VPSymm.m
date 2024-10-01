close all
clear all
clc

rng(936)%rng(693)
v = 23;% v23 good with UB not UB1
suffix = ['fit2VPSymm_multistrat_V' num2str(v)]; %
fit = 0; % 0 if you wanto just to check parameters and simulate
check = 0;
simulate = 1;
simulate_groups = 0;
mode ='dyad';%'single';%'
gen_action_mode = 'stochastic';%'deterministic';%
n_sim = 1; % number of simulations per dyad

%% Getting Dataset
cd ../..
rootfolder = cd;
datadir = fullfile(rootfolder,'data','experiments','expVerticalVPs');
cd src/2VPgame/

nnodes_dataset = 5;

%% Exp Design
n_groups = 2;
n_dyads_group = 9;
n_dyads = n_groups*n_dyads_group;
subH = 1:9; subVH = 10:18;

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
        [tcE,tcL, duration, yrange, start, final, VP1, VP2, tgterror,VPerror] = get_2VPtask_specs('Symm');

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
            x0(n*2-1) =  xmin_jeark(n);
            x0_LB((n*2-1):n*2) = [xmin_jeark(n)-0.0001 -0.0001];%[xmin_jeark(n)-0.01 -0.01];%yrange(1)/2];%[0.1*x01(n*2-1) -0.05];
            x0_UB((n*2-1):n*2) = [xmin_jeark(n)+0.0001 0.0001];%[xmin_jeark(n)+0.01 0.01];%yrange(2)/2];%[10*x01(n*2-1) 0.05];
        end

        % v12
        Xtrue = [0.99,5e-6,1e-2,x0',5e-6,1,0.99];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,5e-5,1e-1,x0_UB,5e-5,2,1-eps];
        
        % v13
        Xtrue = [0.99,5e-6,1e-2,x0',5e-6,1,0.99];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,1e-5,5e-2,x0_UB,1e-5,2,1-eps];

        % v14
        Xtrue = [0.99,5e-6,1e-2,x0',5e-6,1,0.99];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,1e-5,1e-1,x0_UB,1e-5,2,1-eps];
        
        % v15
        Xtrue = [0.99,1e-6,1e-2,x0',1e-6,.1,0.99];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,2e-5,5e-2,x0_UB,2e-5,2,1-eps];
        
        % v16
        Xtrue = [0.999,2e-6,3e-2,x0',2e-6,.1,0.99];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,5e-5,5e-2,x0_UB,5e-5,1,1-eps];
       
        % v17
        Xtrue = [0.999,2e-6,3e-2,x0',2e-6,.1,0.999];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,2e-5,1e-1,x0_UB,1e-5,1,1-eps];
       
        % v18
        Xtrue = [0.999,2e-6,3e-2,x0',2e-6,.1,0.999];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,2e-5,1e-1,x0_UB,2e-5,2,1-eps];

        % v19 nope
        Xtrue = [0.999,2e-6,3e-2,x0',2e-6,.1,0.999];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,1e-4,1e-1,x0_UB,1e-4,2,1-eps];
       
        % v20 
        Xtrue = [0.999,2e-6,3e-2,x0',2e-6,.1,0.999];
        LB = [0.8,eps,eps,x0_LB,eps,eps,0.8];
        UB = [1-eps,1e-5,1e-2,x0_UB,1e-5,0.5,1-eps];

        % v21 improving
        Xtrue = [0.999,2e-6,3e-2,x0',2e-6,1,0.999];
        LB = [0.98,eps,eps,x0_LB,eps,eps,0.98];
        UB = [1-eps,5e-5,1e-1,x0_UB,1e-5,2.5,1-eps];

        % H the other is doing my task
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
        x01_LB = x01'-5e-3;
        x01_UB = x01'+5e-3;

        x02_LB = x02'-5e-3;
        x02_UB = x02'+5e-3;

        % v22 nope
        Xtrue1 = [0.99,2e-6,3e-2,x01',2e-6,1,0.99];
        LB1 = [0.9,eps,eps,x01_LB,eps,eps,0.9];
        UB1 = [1-eps,5e-5,1e-1,x01_UB,1e-5,2.5,1-eps];

        Xtrue2 = [0.99,2e-6,1e-2,x02',2e-6,1.5,0.99];
        LB2 = [0.9,eps,eps,x02_LB,eps,eps,0.9];
        UB2 = [1-eps,5e-5,1e-1,x02_UB,1e-5,2.5,1-eps];

        % % v23 not bad params also sigmax sign with opposite trende wrt
        % % sigmay
        Xtrue1 = [0.999,2e-5,3e-2,x01',2e-6,2,0.999];
        LB1 = [0.9,eps,eps,x01_LB,eps,eps,0.9];
        UB1 = [1-eps,5e-5,1e-1,x01_UB,1e-5,2.5,1-eps];

        Xtrue2 = [0.999,2e-5,1e-2,x02',2e-6,2,0.999];
        LB2 = [0.9,eps,eps,x02_LB,eps,eps,0.9];
        UB2 = [1-eps,5e-5,1e-1,x02_UB,1e-5,2.5,1-eps];
        
        % V24 as V49 Asymm not fitting well
        Xtrue1 = [0.999,1e-6,5e-2,x01',1e-6,.1,0.999];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,2e-6,5e-2,x01_UB,2e-5,1.5,1-eps];
        
        Xtrue2 =  [0.999,1e-6,5e-2,x02',1e-6,.1,0.99];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,2e-6,5e-2,x02_UB,2e-5,1.5,1-eps];

        % v25 nope
        Xtrue1 = [0.999,2e-5,3e-2,x01',2e-6,2,0.999];
        LB1 = [0.9,eps,eps,x01_LB,eps,eps,0.9];
        UB1 = [1-eps,5e-2,1e-2,x01_UB,1e-3,2.5,1-eps];

        Xtrue2 = [0.999,2e-5,1e-2,x02',2e-6,2,0.999];
        LB2 = [0.9,eps,eps,x02_LB,eps,eps,0.9];
        UB2 = [1-eps,5e-2,1e-2,x02_UB,1e-3,2.5,1-eps];
       
        % % v26 as v23 but more memory 
        Xtrue1 = [0.99,2e-5,3e-2,x01',2e-6,2,0.99];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,5e-5,1e-1,x01_UB,1e-5,2.5,1-eps];

        Xtrue2 = [0.99,2e-5,1e-2,x02',2e-6,2,0.99];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,5e-5,1e-1,x02_UB,1e-5,2.5,1-eps];
        
        % % v28 as v26 but more memory and sigmax
        Xtrue1 = [0.98,5e-5,1e-2,x01',2e-6,2,0.98];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,1e-4,1e-1,x01_UB,1e-5,2.5,1-eps];

        Xtrue2 = [0.98,5e-5,1e-2,x02',2e-6,2,0.98];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,1e-4,1e-1,x02_UB,1e-5,2.5,1-eps];

        if d == 10
            Xtrue1 = [0.99,5e-5,1e-2,x01',2e-6,2,0.99];
            Xtrue2 = [0.99,5e-5,1e-2,x02',2e-6,2,0.99];
        end

        % v29 as v28 but less sigmax
        Xtrue1 =[0.98,1e-5,1e-2,x01',2e-6,2,0.98];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,5e-5,1e-1,x01_UB,1e-5,2.5,1-eps];

        Xtrue2 = [0.98,1e-5,1e-2,x02',2e-6,2,0.98];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,5e-5,1e-1,x02_UB,1e-5,2.5,1-eps];

        if d == 10
            Xtrue1 = [0.98,1e-5,1e-2,x01',2e-6,2,0.98];
            Xtrue2 = [0.98,1e-5,1e-2,x02',2e-6,2,0.98];
        end

        % V30 = v28 different init and x0sigmax
        Xtrue1 = [0.9,5e-5,1e-2,x01',2e-6,1,0.98];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,1e-4,1e-1,x01_UB,1e-5,1.5,1-eps];

        Xtrue2 = [0.9,5e-5,1e-2,x02',2e-6,1,0.98];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,1e-4,1e-1,x02_UB,1e-5,1.5,1-eps];


        % V31 = v30 different p0 range
        Xtrue1 = [0.9,5e-5,1e-2,x01',2e-6,1,0.98];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,1e-4,1e-1,x01_UB,1e-4,1.5,1-eps];

        Xtrue2 = [0.9,5e-5,1e-2,x02',2e-6,1,0.98];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,1e-4,1e-1,x02_UB,1e-4,1.5,1-eps];

        % V32 = v31 greater lambda range
        Xtrue1 = [0.9,5e-5,1e-2,x01',2e-6,1,0.98];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,1e-4,1e-1,x01_UB,1e-4,3,1-eps];

        Xtrue2 = [0.9,5e-5,1e-2,x02',2e-6,1,0.98];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,1e-4,1e-1,x02_UB,1e-4,3,1-eps];

        % V33 = v30 smaller a range a
        Xtrue1 = [0.9,5e-5,1e-2,x01',2e-6,1,0.999];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.99];
        UB1 = [1-eps,1e-4,1e-1,x01_UB,1e-5,1.5,1-eps];

        Xtrue2 = [0.9,5e-5,1e-2,x02',2e-6,1,0.999];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.99];
        UB2 = [1-eps,1e-4,1e-1,x02_UB,1e-5,1.5,1-eps];

        % V34 = v30 smaller A range 
        Xtrue1 = [0.999,5e-5,1e-2,x01',2e-6,1,0.99];
        LB1 = [0.99,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,1e-4,1e-1,x01_UB,1e-5,1.5,1-eps];

        Xtrue2 = [0.9,5e-5,1e-2,x02',2e-6,1,0.99];
        LB2 = [0.99,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,1e-4,1e-1,x02_UB,1e-5,1.5,1-eps];

        % V35 = v30 different init and range A
        Xtrue1 = [0.85,5e-5,1e-2,x01',2e-6,1,0.98];
        LB1 = [0.6,eps,eps,x01_LB,eps,eps,0.8];
        UB1 = [1-eps,1e-4,1e-1,x01_UB,1e-5,1.5,1-eps];

        Xtrue2 = [0.85,5e-5,1e-2,x02',2e-6,1,0.98];
        LB2 = [0.6,eps,eps,x02_LB,eps,eps,0.8];
        UB2 = [1-eps,1e-4,1e-1,x02_UB,1e-5,1.5,1-eps];

        % v36 like v23 but more sigma x
        Xtrue1 = [0.999,2e-5,3e-2,x01',2e-6,2,0.999];
        LB1 = [0.9,eps,eps,x01_LB,eps,eps,0.9];
        UB1 = [1-eps,1e-4,1e-1,x01_UB,1e-5,2.5,1-eps];

        Xtrue2 = [0.999,2e-5,1e-2,x02',2e-6,2,0.999];
        LB2 = [0.9,eps,eps,x02_LB,eps,eps,0.9];
        UB2 = [1-eps,1e-4,1e-1,x02_UB,1e-5,2.5,1-eps];
        
        % v37 like v36
        Xtrue1 = [0.999,2e-5,3e-2,x01',2e-6,2,0.999];
        LB1 = [0.8,eps,eps,x01_LB,eps,eps,0.9];
        UB1 = [1-eps,1e-4,1e-1,x01_UB,1e-5,2.5,1-eps];

        Xtrue2 = [0.999,2e-5,1e-2,x02',2e-6,2,0.999];
        LB2 = [0.8,eps,eps,x02_LB,eps,eps,0.9];
        UB2 = [1-eps,1e-4,1e-1,x02_UB,1e-5,2.5,1-eps];
        
        % v38 widen sigmax
        Xtrue1 = [0.999,2e-5,3e-2,x01',2e-6,2,0.999];
        LB1 = [0.6,eps,eps,x01_LB,eps,eps,0.9];
        UB1 = [1-eps,5e-4,1e-1,x01_UB,1e-5,2.5,1-eps];

        Xtrue2 = [0.999,2e-5,1e-2,x02',2e-6,2,0.999];
        LB2 = [0.6,eps,eps,x02_LB,eps,eps,0.9];
        UB2 = [1-eps,5e-4,1e-1,x02_UB,1e-5,2.5,1-eps];
        

        % % v39 like v23 more memory
        Xtrue1 = [0.999,2e-5,3e-2,x01',2e-6,2,0.999];
        LB1 = [0.7,eps,eps,x01_LB,eps,eps,0.7];
        UB1 = [1-eps,5e-5,1e-1,x01_UB,1e-5,2.5,1-eps];

        Xtrue2 = [0.999,2e-5,1e-2,x02',2e-6,2,0.999];
        LB2 = [0.7,eps,eps,x02_LB,eps,eps,0.7];
        UB2 = [1-eps,5e-5,1e-1,x02_UB,1e-5,2.5,1-eps];
        

        X01 = Xtrue1;
        X02 = Xtrue2;

        opts = optimoptions('fmincon');
        %opts.Display = 'iter';

        % partner 1
        X1 = fmincon(@(x) costfunction(x,task1,U2',U1'),X01,[],[],[],[],LB1,UB1,[],opts);
        % partner 2
        X2 = fmincon(@(x) costfunction(x,task2,U1',U2'),X02,[],[],[],[],LB2,UB2,[],opts);

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

        [corrpl1,corrpl2] = get_corr2(nodes_num,U1',uhat1_t,U2',uhat2_t);

        % Saving
        save(fullfile(datadir, ['Dyads_' num2str(d) '_' suffix '.mat']),...
            "Pl1_Params","Pl2_Params","parnames_saving","tcE","tcL","cost1","cost2","U1hat","Phat_U1","pi1","U2hat","Phat_U2","pi2",...
            "dist1","dist2","corrpl1","corrpl2");
    end
    pause
end
%% Checking Parameters After fitting

if check
    GROUPcol = [0.7 0.3 0.7; 0.3 0.7 0.7 ];
    sz = 30;

    % Loading Parameters
    for d = 1:n_dyads
        data_param = load(fullfile(datadir,['Dyads_' num2str(d) '_' suffix '.mat']));
        for p = 1:length(data_param.parnames_saving)
            param = data_param.parnames_saving{p}
            eval([param '1(d,:) = data_param.Pl1_Params.' param ';']);
            eval([param '2(d,:) = data_param.Pl2_Params.' param ';']);
            if p == 4
            else
                eval(['MATRIX_PARAMS((d-1)*2+1,p) = ' param  '1(d,:);']);
                eval(['MATRIX_PARAMS(d*2,p) =' param '2(d,:);']);
                MATRIX_PL([(d-1)*2+1 d*2]) = [1 2];
                MATRIX_DD([(d-1)*2+1 d*2]) = [d d];
            end
        end

        

        nodes_num = size(data_param.U1hat{1,1},1)/2;
        [tcE,tcL, duration, yrange, start, final, VP1, VP2, tgterror,VPerror] = get_2VPtask_specs('Symm');
        twoVPg = twoVP_game(nodes_num,VP1,VP2,tcE,tcL,start,final,duration,VPerror,tgterror);

        for t = 1:size(data_param.pi1,1)
            uhat1_t(:,t) = zeros(nodes_num*2,1);
            uhat2_t(:,t) = zeros(nodes_num*2,1);

            for s = 1:size(data_param.pi1,2)
                uhat1_t(:,t) = uhat1_t(:,t) + data_param.pi1(t,s).*data_param.U1hat{t,s};
                uhat2_t(:,t) = uhat2_t(:,t) + data_param.pi2(t,s).*data_param.U2hat{t,s};
            end
        end

       % plot_minimum_distance(twoVPg,uhat1_t,uhat2_t);
        % [md11_fit(d,:),md12_fit(d,:),md21_fit(d,:),md22_fit(d,:),...
        %     tc11_fit(d,:),tc12_fit(d,:),tc21_fit(d,:),tc22_fit(d,:)] = get_minimum_distance(twoVPg,uhat1_t,uhat2_t);


    end

    MATRIX_PARAMS(:,4) = [];

    % BoxPlot for each parameter
    gr_labels =  {'H','H','H','H','H','H','H','H','H',...
        'VH','VH','VH','VH','VH','VH','VH','VH','VH',...
        'H','H','H','H','H','H','H','H','H',...
        'VH','VH','VH','VH','VH','VH','VH','VH','VH',...
        };
    params = {'A','SigmaX','SigmaY','lambda0','a','P0'};
    tit_params = {'A','$\Sigma x$','$\Sigma y$','$\lambda$','$a$','$\hat{P_0}$'};

    H = [1:9 19:27];
    VH = [10:18 28:36];

    for pp = 1:length(params)

        eval(['p1 =' params{pp} num2str(1) ';']);
        eval(['p2 =' params{pp} num2str(2) ';']);

        eval(['y = [' params{pp} '1; ' params{pp} '2];'])

        % stat
        isnormgr1 = kstest(y(H));
        isnormgr2 = kstest(y(VH));
        if isnormgr1 && isnormgr2
            [p,tbl,stats] = anova1(y,gr_labels);
            [h,p,tbl,stats] = ttest2(y(H),y(VH));
        else
            [p,tbl,stats] = ranksum(y(H),y(VH));
            [p,tbl,stats] = signrank(y(H),y(VH));
        end
        title(params{pp})
        xticklabels({'H','VH'})
        set(gca,'fontname','Times')
        %[results,m,h] = multcompare(stats);
        %tbl = array2table(results,"VariableNames", ...
        %    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


        % boxplot
        maxval = max([p1; p2]);
        minval = min([p1; p2]);

        figure
        set(gcf,'pos',[0 0 125 210])
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
      
        plot(ones(1,18)+(rand(size(ones(1,18)))-0.5).*0.5,[p1(1:9)' p2(1:9)'],'.','color',GROUPcol(2,:),'markersize',10)
        plot(2.*ones(1,18)+(rand(size(ones(1,18)))-0.5).*0.5,[p1(10:18)' p2(10:18)'],'.','color',GROUPcol(1,:),'markersize',10)


        box off
        xticklabels({'H','VH'})
        %title(['p = ' num2str(p)])
        ylabel(tit_params{pp},'interpreter','latex')
        ylim([minval-0.01*minval maxval+0.01*minval])
        set(gca,'fontname','Times')
    end

    % Figure SigmaY-Lambda
    % figure
    % hold on
    % scatter([SigmaY1(subH)' SigmaY2(subH)'],[lambda01(subH)' lambda02(subH)'],sz,GROUPcol(2,:),'filled')
    % scatter([SigmaY1(subVH)' SigmaY2(subVH)'],[lambda01(subVH)' lambda02(subVH)'],sz,GROUPcol(1,:),'filled')
    % ylabel('Lambda')
    % xlabel('Sigma_y')
    % 
    % set(gcf,'pos',[0 0 300 300])


    % Analysis based on learners/non learners
    load(fullfile(datadir,'learners1.mat'))
    load(fullfile(datadir,'learners2.mat'))
    load(fullfile(datadir,'nolearners1.mat'))
    load(fullfile(datadir,'nolearners2.mat'))

    learners = [learners1; learners2+ length(gr_labels)/2];
    nolearners = [nolearners1; nolearners2+ length(gr_labels)/2];
    for ii = 1:length(learners)
        gr_labels{learners(ii)} = 'Learner';
    end
    for ii = 1:length(nolearners)
        gr_labels{nolearners(ii)} = 'NonLearner';
    end

    LEARNcol = [59 179 0]./255;%[0 1 0];
    NOLEARNcol = [255 102 102]./255;%[1 0.5 0];
    GROUPcol = [LEARNcol;NOLEARNcol];

    for pp = 1:length(params)
        eval(['p1 =' params{pp} num2str(1) ';']);
        eval(['p2 =' params{pp} num2str(2) ';']);

        params{pp}
        eval(['y = [' params{pp} '1; ' params{pp} '2];'])

        eval([params{pp} 'av_1 = median(y(learners));' ])
        eval([params{pp} 'perc40_1 = prctile(y(learners),45);' ])
        eval([params{pp} 'perc60_1 = prctile(y(learners),55);' ])
        eval([params{pp} 'av_2 = median(y(nolearners));' ])
        eval([params{pp} 'perc40_2 = prctile(y(nolearners),45);' ])
        eval([params{pp} 'perc60_2 = prctile(y(nolearners),55);' ])
        eval([params{pp} 'std_1 = std(y(learners));' ])
        eval([params{pp} 'std_2 = std(y(nolearners));' ])

        isnormgr1 = kstest(y(learners));
        isnormgr2 = kstest(y(nolearners));
        if isnormgr1 && isnormgr2
            [p,tbl,stats] = anova1(y,gr_labels);
            [h,p,tbl,stats] = ttest2(y(learners),y(nolearners)); % not balanced
        else
            [p,tbl,stats] = ranksum(y(learners),y(nolearners));
            [p,tbl,stats] = signrank(y(learners),y(nolearners));
        end
        title(params{pp})
        xticklabels({'Learn','No-Learn'})
        set(gca,'fontname','Times')

        maxval = max([p1; p2]);
        minval = min([p1 ;p2]);

        figure
        set(gcf,'pos',[0 0 150 250])
        % boxplot([reshape(p1,n_dyads,n_groups);reshape(p2,n_dyads,n_groups)])
        boxplot(y,gr_labels)

        h = findobj(gca,'Tag','Box');
        mm = findobj(gca,'Tag','Median');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),GROUPcol(j,:),'FaceAlpha',.3,...
                'Edgecolor',GROUPcol(j,:));
        end

        hold on
        gr_bool(learners) = 1; gr_bool(nolearners) = 0;

        plot(ones(1,length(learners))+(rand(size(ones(1,length(learners))))-0.5).*0.5,[y(learners)],'.','color',GROUPcol(2,:),'markersize',10)
        plot(2.*ones(1,length(nolearners))+(rand(size(ones(1,length(nolearners))))-0.5).*0.5,[y(nolearners)],'.','color',GROUPcol(1,:),'markersize',10)
        box off
        xticklabels({'Learners','Non Learners'})
        %title(['p = ' num2str(p)])
        ylabel(tit_params{pp},'Interpreter','latex')
        ylim([minval-0.01*minval maxval+0.01*minval])
        set(gca,'fontname','Times')
    end

    %GROUPcol = [0.7 0.3 0.7; 0.3 0.7 0.7 ];
    params = {'A','SigmaX','SigmaY','lambda0','a','P0'};
    tit_params = {'A','$\Sigma x$','$\Sigma y$','$\lambda$','$a$','$\hat{P_0}$'};
    for pp = 1:length(params)
        params{pp}

        eval(['p1 =' params{pp} num2str(1) ';']);
        eval(['p2 =' params{pp} num2str(2) ';']);

        learnersH = [p1(intersect(subH,learners1)); p2(intersect(subH,learners2))];
        learnersVH = [p1(intersect(subVH,learners1)); p2(intersect(subVH,learners2))];

        nonlearnersH = [p1(intersect(subH,nolearners1)); p2(intersect(subH,nolearners2))];
        nonlearnersVH = [p1(intersect(subVH,nolearners1)); p2(intersect(subVH,nolearners2))];

        %figname = [tit_params{pp} '.fig'];
        xgroupdata = {ones(size(learnersH)), 3.5.*ones(size(nonlearnersH));...
            2.*ones(size(learnersVH)), 4.5.*ones(size(nonlearnersVH))};
        ydata = {learnersH nonlearnersH ; learnersVH  nonlearnersVH;};
        colorgroup = { GROUPcol(2,:) GROUPcol(1,:);...
            GROUPcol(2,:) GROUPcol(1,:)};
        col =colorgroup;
        xticklab = {'H','VH','H','VH'};
        groupbox_model('ff.fig',xgroupdata,ydata,colorgroup,col,xticklab,tit_params{pp});
        ylabel(tit_params{pp},'Interpreter','Latex');
        %title(tit,'Interpreter','Latex')
        set(gca,'fontname','times')
        xlim([0.5 5])

        set(gcf,'Position',[0 0 200 250])
    end

    % % PLot Fitting Minimum Distance
    % colorgroup = [0.3 0.7 0.7; 0.7 0.3 0.7; 0.7 0.7 0.3 ];
    % alpha = .3;
    % gg = {subH,subVH};
    % n_sub = n_dyads_group;
    % T = size(md11_fit,2);
    % %colorgroup = GROUPcol;
    % 
    % 
    % 
    % mdlearnersH1 = md12_fit(intersect(subH,learners1),:);
    % n_sublH1 = size(mdlearnersH1,1);
    % mdnolearnersH1 = md12_fit(intersect(subH,nolearners1),:);
    % n_subnlH1 = size(mdnolearnersH1,1);
    % 
    % mdlearnersVH1 = md12_fit(intersect(subVH,learners1),:);
    % n_sublVH1 = size(mdlearnersVH1,1);
    % mdnolearnersVH1 = md12_fit(intersect(subVH,nolearners1),:);
    % n_subnlVH1 = size(mdnolearnersVH1,1);
    % 
    % mdlearnersH2 = md21_fit(intersect(subH,learners2),:);
    % n_sublH2 = size(mdlearnersH2,1);
    % mdnolearnersH2 = md21_fit(intersect(subH,nolearners2),:);
    % n_subnlH2 = size(mdnolearnersH2,1);
    % 
    % mdlearnersVH2 = md21_fit(intersect(subVH,learners2),:);
    % n_sublVH2 = size(mdlearnersVH2,1);
    % mdnolearnersVH2 = md21_fit(intersect(subVH,nolearners2),:);
    % n_subnlVH2 = size(mdnolearnersVH2,1);
    % 
    % figure
    % subplot(121)
    % hold on
    % patch([1:T T:-1:1],[mean(mdlearnersH1,1)-std(mdlearnersH1,1)./sqrt(n_sublH1) fliplr(mean(mdlearnersH1,1)+std(mdlearnersH1,1)./sqrt(n_sublH1))],...
    %     colorgroup(1,:),'facealpha',alpha,'edgecolor',colorgroup(1,:),'edgealpha',alpha)
    % 
    % patch([1:T T:-1:1],[mean(mdlearnersVH1,1)-std(mdlearnersVH1,1)./sqrt(n_subnlH1) fliplr(mean(mdlearnersVH1,1)+std(mdlearnersVH1,1)./sqrt(n_subnlVH1))],...
    %     colorgroup(2,:),'facealpha',alpha,'edgecolor',colorgroup(2,:),'edgealpha',alpha)
    % 
    % s(1) = plot(mean(mdlearnersH1,1),'linewidth',1.5,'color',colorgroup(1,:))
    % s(2) = plot(mean(mdlearnersVH1,1),'linewidth',1.5,'color',colorgroup(2,:))
    % 
    % ylim([0 0.06])
    % xlim([0 T+1])
    % xlabel('Trial')
    % ylabel('distance [m]')
    % 
    % legend(s,{'H','VH'})
    % legend boxoff
    % title('MD12 Learners')
    % 
    % set(gca,'fontname','times')
    % 
    % subplot(122)
    % hold on
    % patch([1:T T:-1:1],[mean(mdnolearnersH1,1)-std(mdnolearnersH1,1)./sqrt(n_subnlH1) fliplr(mean(mdnolearnersH1,1)+std(mdnolearnersH1,1)./sqrt(n_subnlH1))],...
    %     colorgroup(1,:),'facealpha',alpha,'edgecolor',colorgroup(1,:),'edgealpha',alpha)
    % 
    % patch([1:T T:-1:1],[mean(mdnolearnersVH1,1)-std(mdnolearnersVH1,1)./sqrt(n_subnlVH1) fliplr(mean(mdnolearnersVH1,1)+std(mdnolearnersVH1,1)./sqrt(n_subnlVH1))],...
    %     colorgroup(2,:),'facealpha',alpha,'edgecolor',colorgroup(2,:),'edgealpha',alpha)
    % 
    % s(1) = plot(mean(mdnolearnersH1,1),'linewidth',1.5,'color',colorgroup(1,:))
    % s(2) = plot(mean(mdnolearnersVH1,1),'linewidth',1.5,'color',colorgroup(2,:))
    % 
    % ylim([0 0.06])
    % xlim([0 T+1])
    % xlabel('Trial')
    % ylabel('distance [m]')
    % 
    % legend(s,{'H','VH'})
    % legend boxoff
    % title('MD12 No Learners')
    % set(gca,'fontname','times')
    % 
    % set(gcf,'pos',[0 0 300 190])
    % 
    % 
    % figure
    % subplot(121)
    % hold on
    % patch([1:T T:-1:1],[mean(mdlearnersH2,1)-std(mdlearnersH2,1)./sqrt(n_sublH2) fliplr(mean(mdlearnersH2,1)+std(mdlearnersH2,1)./sqrt(n_sublH2))],...
    %     colorgroup(1,:),'facealpha',alpha,'edgecolor',colorgroup(1,:),'edgealpha',alpha)
    % 
    % patch([1:T T:-1:1],[mean(mdlearnersVH2,1)-std(mdlearnersVH2,1)./sqrt(n_subnlH2) fliplr(mean(mdlearnersVH2,1)+std(mdlearnersVH2,1)./sqrt(n_subnlVH2))],...
    %     colorgroup(2,:),'facealpha',alpha,'edgecolor',colorgroup(2,:),'edgealpha',alpha)
    % 
    % s(1) = plot(mean(mdlearnersH2,1),'linewidth',1.5,'color',colorgroup(1,:))
    % s(2) = plot(mean(mdlearnersVH2,1),'linewidth',1.5,'color',colorgroup(2,:))
    % 
    % ylim([0 0.06])
    % xlim([0 T+1])
    % xlabel('Trial')
    % ylabel('distance [m]')
    % 
    % legend(s,{'H','VH'})
    % legend boxoff
    % title('MD21 Learners')
    % 
    % set(gca,'fontname','times')
    % 
    % subplot(122)
    % hold on
    % patch([1:T T:-1:1],[mean(mdnolearnersH2,1)-std(mdnolearnersH2,1)./sqrt(n_subnlH2) fliplr(mean(mdnolearnersH2,1)+std(mdnolearnersH2,1)./sqrt(n_subnlH2))],...
    %     colorgroup(1,:),'facealpha',alpha,'edgecolor',colorgroup(1,:),'edgealpha',alpha)
    % 
    % patch([1:T T:-1:1],[mean(mdnolearnersVH2,1)-std(mdnolearnersVH2,1)./sqrt(n_subnlVH1) fliplr(mean(mdnolearnersVH2,1)+std(mdnolearnersVH2,1)./sqrt(n_subnlVH2))],...
    %     colorgroup(2,:),'facealpha',alpha,'edgecolor',colorgroup(2,:),'edgealpha',alpha)
    % 
    % s(1) = plot(mean(mdnolearnersH2,1),'linewidth',1.5,'color',colorgroup(1,:))
    % s(2) = plot(mean(mdnolearnersVH2,1),'linewidth',1.5,'color',colorgroup(2,:))
    % 
    % ylim([0 0.06])
    % xlim([0 T+1])
    % xlabel('Trial')
    % ylabel('distance [m]')
    % 
    % legend(s,{'H','VH'})
    % legend boxoff
    % title('MD21 No Learners')
    % set(gca,'fontname','times')
    % 
    % set(gcf,'pos',[0 0 300 190])
    % 

MATRIX = [];

MATRIX = [(1:36)' MATRIX_DD' MATRIX_PL' MATRIX_PARAMS];
save(fullfile('/Users/ceciliadevicariis/Documents/MATLAB/dati_ViaPoints_Shift_Symm/SYM_INDEXES',['MATRIX_PARAMSV' num2str(v) '.mat']),"MATRIX");
pause

end



%% Simulation dyad by dyad
if simulate

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
        [tcE,tcL, duration, yrange, start, final, VP1, VP2, tgterror,VPerror] = get_2VPtask_specs('Symm');

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

            % plot trajectory and partner model animation
            %plot_trajectory(twoVPg,u1,u2,u1ne,u2ne)
            %plot_partnermodel_action(twoVPg,x1_prior,x2_prior,u1,u2)

            % plot actions
            %plot_action(twoVPg,u1(:,1:10),u2(:,1:10));
            figure
            subplot(121)
            plot_action(twoVPg,U1((T-9):T,:)',U2((T-9):T,:)');

            subplot(122)
            plot_action(twoVPg,u1(:,(end-9):end),u2(:,(end-9):end));
            % figure
            % plot_action(twoVPg,u1,u2);
            set(gcf,'pos',[0 500 250 125])
            % if s == 1
            %     figure
            %     plot_action(twoVPg,U1((T-10):T,:)',U2((T-10):T,:)');
            %     set(gcf,'pos',[250 500 250 250])
            % 
            % end

            %plot_animation(twoVPg,u1,u2);

            % plot speed profile
            %plot_speed(twoVPg,u1,u2);

            % plot acceleration profile
            %plot_acceleration(twoVPg,u1,u2);

            % plot minimum distance
            % plot_minimum_distance(twoVPg,u1,u2);
            [md11_sim{d}(s,:),md12_sim{d}(s,:),md21_sim{d}(s,:),md22_sim{d}(s,:),...
                tc11_sim{d}(s,:),tc12_sim{d}(s,:),tc21_sim{d}(s,:),tc22_sim{d}(s,:)] = get_minimum_distance(twoVPg,u1,u2);

            if s == 1
                exp_sig = [md11(d,:);md12(d,:);md21(d,:);md22(d,:)];
                sim_sig = [md11_sim{d}(s,:);md12_sim{d}(s,:);md21_sim{d}(s,:);md22_sim{d}(s,:)];
                % plot_md_exp_sim(twoVPg,exp_sig,sim_sig,1:T,[0 0.05],{'MD_{11}','MD_{12}','MD_{21}','MD{22}'})
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
            % plot_traj_differences(twoVPg,nodes_num,U1',u1,U2',u2);

            [corrpl1{d}(s,:),corrpl2{d}(s,:)] = get_corr2(nodes_num,U1',u1,U2',u2);

            %% Symm Collaborative behavior classification
            th_md = 0.02;
            ACTION_pl1{d}(s,:) = get_action_class(th_md,md12_sim{d}(s,:),tc11_sim{d}(s,:),tc12_sim{d}(s,:));
            ACTION_pl2{d}(s,:) = get_action_class(th_md,md21_sim{d}(s,:),tc21_sim{d}(s,:),tc22_sim{d}(s,:));

            %% Probabilities
            [collab{d}(s,:),collab12{d}(s,:),collab21{d}(s,:),oscil{d}(s,:),nocollab{d}(s,:),oscil12_21{d}(s,:),oscil21_12{d}(s,:),bothignore{d}(s,:),oneignore{d}(s,:)] = get_prob_ja(ACTION_pl1{d}(s,:),ACTION_pl2{d}(s,:))

            u1temp = u1sim{d};
            u2temp = u2sim{d};
            dist1temp = dist1_sim{d}; dist2temp = dist2_sim{d};
            spvar1temp = spvar1{d}; spvar2temp = spvar2{d};
            corrcoeff1temp = corrpl1{d}; corrcoeff2temp = corrpl2{d};
            save(fullfile(rootfolder,'data','simulations','VPSymm',['dyad_' num2str(d) '_V' num2str(v) '_' mode '.mat']),...
                'U1','U2','u1temp','u2temp',"dist1temp","dist2temp", "spvar2temp","spvar2temp",...
                "corrcoeff1temp",'corrcoeff2temp')

        end
        close all
    end


    % Plot Group MDs H-VH
    colorgroup = [0.3 0.7 0.7; 0.7 0.3 0.7; 0.7 0.7 0.3 ];
    alpha = .3
    indicators = {};%'md12_sim','md21_sim'}
    gg = {subH,subVH};
   
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
        end

    end
    set(gcf,'pos',[0 0 400 210])

    % Plot Group MDs Learners non learners
    load('learners1.mat')
    load('learners2.mat')
    load('nolearners1.mat')
    load('nolearners2.mat')
    LEARNcol = [59 179 0]./255;%[0 1 0];
    NOLEARNcol = [255 102 102]./255;%[1 0.5 0];
    colorgroup = [NOLEARNcol; LEARNcol; 0.7 0.7 0.3 ];
    alpha = .3

    if n_sim>1
        for d = 1:n_dyads
            md12_sim{d} = mean(md12_sim{d});
            md21_sim{d} = mean(md21_sim{d});
        end
    end
    md12_sim_new = cell2mat(md12_sim);
    md21_sim_new = cell2mat(md21_sim);
    
    md12_sim_new = reshape(md12_sim_new,170,18)';
    md21_sim_new = reshape(md21_sim_new,170,18)';

    mdlearners = [md12_sim_new(learners1,:);md21_sim_new(learners2,:)];
    mdnonlearners = [md12_sim_new(nolearners1,:);md21_sim_new(nolearners2,:)];

    nlearners = size(mdlearners,1);
    nnolearners = size(mdnonlearners,1);
    
    for d = 1:size(mdlearners,1)
        mdlearners(d,:) = movmean(mdlearners(d,:),1);
    end
    for d = 1:size(mdnonlearners,1)
        mdnonlearners(d,:) = movmean(mdnonlearners(d,:),1);
    end

    figure
    hold on
    patch([1:T T:-1:1],[mean(mdlearners,1)-std(mdlearners,1)./sqrt(nlearners) fliplr(mean(mdlearners,1)+std(mdlearners,1)./sqrt(nlearners))],...
        colorgroup(1,:),'facealpha',alpha,'edgecolor',colorgroup(1,:),'edgealpha',alpha)

    patch([1:T T:-1:1],[mean(mdnonlearners,1)-std(mdnonlearners,1)./sqrt(nnolearners) fliplr(mean(mdnonlearners,1)+std(mdnonlearners,1)./sqrt(nnolearners))],...
        colorgroup(2,:),'facealpha',alpha,'edgecolor',colorgroup(2,:),'edgealpha',alpha)

    s(1) = plot(mean(mdlearners,1),'linewidth',1.5,'color',colorgroup(1,:))
    s(2) = plot(mean(mdnonlearners,1),'linewidth',1.5,'color',colorgroup(2,:))

    ylim([0 0.05])
    xlim([0 T+1])
    xlabel('Trial','fontname','time')
    ylabel('distance [m]','fontname','time')
    % title(titles{i})

    legend(s,{'Learners','Non Learners'})
    legend boxoff
    set(gca, 'FontName', 'times');

    set(gcf,'pos',[300 0 150 190])


    % p(a=1) collab; p(a=2) collab; p(a=0)
    purple = [138,43,226]./255;
    orange = [230, 149, 0]./255;
    light_green = [0.1 0.9 0.1];
    green = [0.1 0.6 0.1];
    light_blue = [51 153 255]./255;
    blue = [0 102 204]./255;
    alpha_shade = 0.2;
    lw = 1.5;
    trials_ep = 10;
    n_ep = T/trials_ep;

    % Actions over time
    %     for d = 1:n_dyads
    %     for s = 1:n_sim
    %         figure
    %         set(gcf,'pos',[0 0 190 110])
    %         hold on
    %         stairs(ACTION_pl1{d}(s,:),'color','b','linewidth',lw)
    %         stairs(ACTION_pl2{d}(s,:)+3,'color','r','linewidth',lw)
    %         yticks([0 1 2 3 4 5])
    %         ylim([-0.1 5.1])
    %         yticklabels({'11' '12' '21' '22' '12' '21'})
    %         title(['Actions D' num2str(d)])
    %         xlabel('Trials')
    %         xlim([0 170])
    %     end
    % end

% Actions at dyad level
for d = 1:n_dyads
    for s = 1:n_sim
        for t = 1:T
            actEE{d}(s,t) = (ACTION_pl1{d}(s,t) == 1) && (ACTION_pl2{d}(s,t) == 1);
            actLL{d}(s,t) = (ACTION_pl1{d}(s,t) == 2) && (ACTION_pl2{d}(s,t) == 2);

            actEL{d}(s,t) = (ACTION_pl1{d}(s,t) == 1) && (ACTION_pl2{d}(s,t) == 2);
            actLE{d}(s,t) = (ACTION_pl1{d}(s,t) == 2) && (ACTION_pl2{d}(s,t) == 1);

            actEM{d}(s,t) = (ACTION_pl1{d}(s,t) == 1) && (ACTION_pl2{d}(s,t) == 0);
            actLM{d}(s,t) = (ACTION_pl1{d}(s,t) == 2) && (ACTION_pl2{d}(s,t) == 0);

            actME{d}(s,t) = (ACTION_pl1{d}(s,t) == 0) && (ACTION_pl2{d}(s,t) == 1);
            actML{d}(s,t) = (ACTION_pl1{d}(s,t) == 0) && (ACTION_pl2{d}(s,t) == 2);

            actMM{d}(s,t) = (ACTION_pl1{d}(s,t) == 0) && (ACTION_pl2{d}(s,t) == 0);

        end
        act_EE_or_LL{d}(s,:) = actEE{d}(s,:) + actLL{d}(s,:); % collab
        act_EL_or_LE{d}(s,:) = actEL{d}(s,:) + actLE{d}(s,:); % oscill
        act_EM_LM_ME_LE{d}(s,:) = actEM{d}(s,:) + actLM{d}(s,:) + actME{d}(s,:) + actML{d}(s,:) + actMM{d}(s,:); % 1collab
        act_nocollab{d}(s,:) = act_EL_or_LE{d}(s,:) + act_EM_LM_ME_LE{d}(s,:);

    end
end

% Probabilities
ntrials_ep = 10;
for d = 1:n_dyads
    for s = 1:n_sim
        for e = 1:n_ep
            currtrials = (1:ntrials_ep) + (e-1)*ntrials_ep ;
            pEEorLL{d}(s,e) = mean(act_EE_or_LL{d}(s,currtrials));
            pELorLE{d}(s,e) = mean(act_EL_or_LE{d}(s,currtrials));
            pEMorLMorMEorLE{d}(s,e) = mean(act_EM_LM_ME_LE{d}(s,currtrials));
            pMM{d}(s,e) = mean(actMM{d}(s,currtrials));
            pnocollab{d}(s,e) = mean(act_nocollab{d}(s,currtrials));

        end

        % figure
        % set(gcf,'pos',[0 0 150 150])
        % hold on
        % stairs(pEEorLL{d}(s,:),'Color',orange,'LineWidth',lw)
        % stairs(pELorLE{d}(s,:),'Color','magenta','LineWidth',lw)
        % stairs(pEMorLMorMEorLE{d}(s,:),'Color','black','LineWidth',lw)
        % %stairs(pMM(d,:),'Color','black','LineWidth',lw)
        % title(['P(diadic actions) D' num2str(s)])
        % xlabel('Epochs')
        % xlim([0 17])
        % ylim([0 1])
        % legend('collab','cycling','not collab')
        % legend boxoff
    end
end

% % Population figure
% for d = 1:n_groups
%     figure
%     set(gcf,'pos',[0 0 150 150])
%     hold on
%     patch([1:n_ep n_ep:-1:1],[mean(pEEorLL{d})-std(pEEorLL{d})./sqrt(n_sim) fliplr(mean(pEEorLL{d})+std(pEEorLL{d})./sqrt(n_sim))],...
%         orange,'FaceAlpha',alpha_shade,'EdgeAlpha',0)
%     patch([1:n_ep n_ep:-1:1],[mean(pELorLE{d})-std(pELorLE{d})./sqrt(n_sim) fliplr(mean(pELorLE{d}) + std(pELorLE{d})./sqrt(n_sim))],...
%         'magenta','Facealpha',alpha_shade,'EdgeAlpha',0)
%     patch([1:n_ep n_ep:-1:1],[mean(pEMorLMorMEorLE{d}) - std(pEMorLMorMEorLE{d})./sqrt(n_sim) fliplr(mean(pEMorLMorMEorLE{d}) + std(pEMorLMorMEorLE{d})./sqrt(n_sim))],...
%         'k','Facealpha',alpha_shade,'EdgeAlpha',0)
% 
%     p(1) = plot(mean(pEEorLL{d}),'Color',orange,'LineWidth',lw)
%     p(2) = plot(mean(pELorLE{d}),'Color','magenta','LineWidth',lw)
%     p(3) = plot(mean(pEMorLMorMEorLE{d}),'Color','black','LineWidth',lw)
% 
%     title(['P(diadic actions)'])
%     xlabel('Epochs')
%     xlim([0 17])
%     ylim([0 1])
%     legend(p,{'collab','cycling','not collab'})
%     legend boxoff
% end

% Population figure

learnersdyads = intersect(learners1,learners2);
ndyads_learners = size(learnersdyads,1);
ndyads_nonlearners = 18-ndyads_learners;
n_sub = [ndyads_learners ndyads_nonlearners];

nonlearnersdyads = 1:18;
nonlearnersdyads(learnersdyads) = [];

if n_sim == 1
    pEEorLL_mat = reshape(cell2mat(pEEorLL),17,18)';
    pELorLE_mat = reshape(cell2mat(pELorLE),17,18)';
    pnocollab_mat = reshape(cell2mat(pnocollab),17,18)';
else
    pEEorLL_mat = reshape(mean(cell2mat(pEEorLL)),17,18)';
    pELorLE_mat = reshape(mean(cell2mat(pELorLE)),17,18)';
    pnocollab_mat = reshape(mean(cell2mat(pnocollab)),17,18)';
end

pEEorLL_new{1} = pEEorLL_mat(learnersdyads,:);
pEEorLL_new{2} = pEEorLL_mat(nonlearnersdyads,:);
pnocollab_new{1} = pnocollab_mat(learnersdyads,:);
pnocollab_new{2} = pnocollab_mat(nonlearnersdyads,:);

tits = {'Learners','Non Learners'}
for g = 1:n_groups
    figure
    set(gcf,'pos',[0 0 150 150])
    hold on
    patch([1:n_ep n_ep:-1:1],[mean(pEEorLL_new{g})-std(pEEorLL_new{g})./sqrt(n_sub(g)) fliplr(mean(pEEorLL_new{g})+std(pEEorLL_new{g})./sqrt(n_sub(g)))],...
        purple,'FaceAlpha',alpha_shade,'EdgeAlpha',0)

    patch([1:n_ep n_ep:-1:1],[mean(pnocollab_new{g}) - std(pnocollab_new{g})./sqrt(n_sub(g)) fliplr(mean(pnocollab_new{g}) + std(pnocollab_new{g})./sqrt(n_sub(g)))],...
        light_green,'Facealpha',alpha_shade,'EdgeAlpha',0)

    p(1) = plot(mean(pEEorLL_new{g}),'Color',purple,'LineWidth',lw)
    p(2) = plot(mean(pnocollab_new{g}),'Color',light_green,'LineWidth',lw)

    title(tits{g})
    xlabel('Epochs')
    xlim([0 17])
    ylim([0 1])
    legend(p,{'Collaborative','Solo'})
    legend boxoff
    legend off
    ylabel('Probaility')
    set(gca,'fontname','times')

end

end




%% Simulation dyad by dyad
if simulate_groups

    clear corrpl1 corrpl2
    load(fullfile(datadir, ['forces.mat']));
    indicators = {'md11','md22','md12','md21','tc11','tc22','te1','te2','ts1','ts2'};
    unconnected = find(forces(end,:) == 0);



    for i = 1:length(indicators)
        load(fullfile(datadir, [indicators{i} '.mat']));
        eval([indicators{i} '= ind_ts;' ]);
        eval([indicators{i} '(:,unconnected)= [];' ]);
    end

    load('learners1')
    load('learners2')
    learners_dyads = intersect(learners1,learners2);
    nonlearners_dyads = 1:18;
    nonlearners_dyads(learners_dyads) = [];
    dyads_grouped = {learners_dyads, nonlearners_dyads};
    Params = defining_player_model_from_fitting_bygroup_symm(datadir,suffix,n_dyads,length(dyads_grouped),dyads_grouped);
    
    for d = 1:n_dyads
        % Defining task specs
        [tcE,tcL, duration, yrange, start, final, VP1, VP2, tgterror,VPerror] = get_2VPtask_specs('Symm');

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

            % plot trajectory and partner model animation
            %plot_trajectory(twoVPg,u1,u2,u1ne,u2ne)
            %plot_partnermodel_action(twoVPg,x1_prior,x2_prior,u1,u2)

            % plot actions
            figure
            subplot(121)
            plot_action(twoVPg,U1((170-10):170,:)',U2((170-10):170,:)');

            subplot(122)
            plot_action(twoVPg,u1(:,(end-9):end),u2(:,(end-9):end));
            % figure
            % plot_action(twoVPg,u1,u2);
            set(gcf,'pos',[0 500 250 125])

            % if s == 1
            %     figure
            %     plot_action(twoVPg,U1((T-10):T,:)',U2((T-10):T,:)');
            %     set(gcf,'pos',[250 500 250 250])
            % 
            % end

            %plot_animation(twoVPg,u1,u2);

            % plot speed profile
            %plot_speed(twoVPg,u1,u2);

            % plot acceleration profile
            %plot_acceleration(twoVPg,u1,u2);

            % plot minimum distance
            % plot_minimum_distance(twoVPg,u1,u2);
            [md11_sim{d}(s,:),md12_sim{d}(s,:),md21_sim{d}(s,:),md22_sim{d}(s,:),...
                tc11_sim{d}(s,:),tc12_sim{d}(s,:),tc21_sim{d}(s,:),tc22_sim{d}(s,:)] = get_minimum_distance(twoVPg,u1,u2);

            if s == 1
                exp_sig = [md11(d,:);md12(d,:);md21(d,:);md22(d,:)];
                sim_sig = [md11_sim{d}(s,:);md12_sim{d}(s,:);md21_sim{d}(s,:);md22_sim{d}(s,:)];
                % plot_md_exp_sim(twoVPg,exp_sig,sim_sig,1:T,[0 0.05],{'MD_{11}','MD_{12}','MD_{21}','MD{22}'})
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
            % plot_traj_differences(twoVPg,nodes_num,U1',u1,U2',u2);

            [corrpl1{d}(s,:),corrpl2{d}(s,:)] = get_corr2(nodes_num,U1',u1,U2',u2);

            %% Symm Collaborative behavior classification
            th_md = 0.02;
            ACTION_pl1{d}(s,:) = get_action_class(th_md,md12_sim{d}(s,:),tc11_sim{d}(s,:),tc12_sim{d}(s,:));
            ACTION_pl2{d}(s,:) = get_action_class(th_md,md21_sim{d}(s,:),tc21_sim{d}(s,:),tc22_sim{d}(s,:));

            %% Probabilities
            [collab{d}(s,:),collab12{d}(s,:),collab21{d}(s,:),oscil{d}(s,:),nocollab{d}(s,:),oscil12_21{d}(s,:),oscil21_12{d}(s,:),bothignore{d}(s,:),oneignore{d}(s,:)] = get_prob_ja(ACTION_pl1{d}(s,:),ACTION_pl2{d}(s,:))

            u1temp = u1sim{d};
            u2temp = u2sim{d};
            dist1temp = dist1_sim{d}; dist2temp = dist2_sim{d};
            spvar1temp = spvar1{d}; spvar2temp = spvar2{d};
            corrcoeff1temp = corrpl1{d}; corrcoeff2temp = corrpl2{d};
            save(fullfile(rootfolder,'data','simulations','VPSymm',['dyad_' num2str(d) '_V' num2str(v) '_' mode '.mat']),...
                'U1','U2','u1temp','u2temp',"dist1temp","dist2temp", "spvar2temp","spvar2temp",...
                "corrcoeff1temp",'corrcoeff2temp')

        end
        %close all
    end


    % Plot Group MDs H-VH
    colorgroup = [0.3 0.7 0.7; 0.7 0.3 0.7; 0.7 0.7 0.3 ];
    alpha = .3
    indicators = {'md12_sim','md21_sim'}
    gg = {subH,subVH};
   
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
        end

    end
    set(gcf,'pos',[0 0 400 210])

    % Plot Group MDs Learners non learners
    load('learners1.mat')
    load('learners2.mat')
    load('nolearners1.mat')
    load('nolearners2.mat')

    LEARNcol = [59 179 0]./255;%[0 1 0];
    NOLEARNcol = [255 102 102]./255;%[1 0.5 0];
    colorgroup = [NOLEARNcol; LEARNcol; 0.7 0.7 0.3 ];
    alpha = .3

    if n_sim>1
        for d = 1:n_dyads
            md12_sim{d} = mean(md12_sim{d});
            md21_sim{d} = mean(md21_sim{d});
        end
    end
    md12_sim_new = cell2mat(md12_sim);
    md21_sim_new = cell2mat(md21_sim);
    
    md12_sim_new = reshape(md12_sim_new,T,18)';
    md21_sim_new = reshape(md21_sim_new,T,18)';

    mdlearners = [md12_sim_new(learners1,:);md21_sim_new(learners2,:)];
    mdnonlearners = [md12_sim_new(nolearners1,:);md21_sim_new(nolearners2,:)];

    nlearners = size(mdlearners,1);
    nnolearners = size(mdnonlearners,1);
    
    for d = 1:size(mdlearners,1)
        mdlearners(d,:) = movmean(mdlearners(d,:),1);
    end
    for d = 1:size(mdnonlearners,1)
        mdnonlearners(d,:) = movmean(mdnonlearners(d,:),1);
    end

    figure
    hold on
    patch([1:T T:-1:1],[mean(mdlearners,1)-std(mdlearners,1)./sqrt(nlearners) fliplr(mean(mdlearners,1)+std(mdlearners,1)./sqrt(nlearners))],...
        colorgroup(1,:),'facealpha',alpha,'edgecolor',colorgroup(1,:),'edgealpha',alpha)

    patch([1:T T:-1:1],[mean(mdnonlearners,1)-std(mdnonlearners,1)./sqrt(nnolearners) fliplr(mean(mdnonlearners,1)+std(mdnonlearners,1)./sqrt(nnolearners))],...
        colorgroup(2,:),'facealpha',alpha,'edgecolor',colorgroup(2,:),'edgealpha',alpha)

    s(1) = plot(mean(mdlearners,1),'linewidth',1.5,'color',colorgroup(1,:))
    s(2) = plot(mean(mdnonlearners,1),'linewidth',1.5,'color',colorgroup(2,:))

    ylim([0 0.05])
    xlim([0 T+1])
    xlabel('Trial','fontname','time')
    ylabel('distance [m]','fontname','time')
    % title(titles{i})

    legend(s,{'Learners','Non Learners'})
    legend boxoff
    set(gca, 'FontName', 'times');

    set(gcf,'pos',[300 0 150 190])

    clear mdlearners mdnonlearners
    % Learners DYADS!
    learnersdyads = intersect(learners1,learners2);
    ndyads_learners = size(learnersdyads,1);
    ndyads_nonlearners = 18-ndyads_learners;
    n_sub = [ndyads_learners ndyads_nonlearners];
    nonlearnersdyads = 1:18;
    nonlearnersdyads(learnersdyads) = [];

    mdlearners = [md12_sim_new(learnersdyads,:);md21_sim_new(learnersdyads,:)];
    mdnonlearners = [md12_sim_new(nonlearnersdyads,:);md21_sim_new(nonlearnersdyads,:)];

    nlearners = size(mdlearners,1);
    nnolearners = size(mdnonlearners,1);
    
    for d = 1:size(mdlearners,1)
        mdlearners(d,:) = movmean(mdlearners(d,:),1);
    end
    for d = 1:size(mdnonlearners,1)
        mdnonlearners(d,:) = movmean(mdnonlearners(d,:),1);
    end

    figure
    hold on
    patch([1:T T:-1:1],[mean(mdlearners,1)-std(mdlearners,1)./sqrt(nlearners) fliplr(mean(mdlearners,1)+std(mdlearners,1)./sqrt(nlearners))],...
        colorgroup(1,:),'facealpha',alpha,'edgecolor',colorgroup(1,:),'edgealpha',alpha)

    patch([1:T T:-1:1],[mean(mdnonlearners,1)-std(mdnonlearners,1)./sqrt(nnolearners) fliplr(mean(mdnonlearners,1)+std(mdnonlearners,1)./sqrt(nnolearners))],...
        colorgroup(2,:),'facealpha',alpha,'edgecolor',colorgroup(2,:),'edgealpha',alpha)

    s(1) = plot(mean(mdlearners,1),'linewidth',1.5,'color',colorgroup(1,:))
    s(2) = plot(mean(mdnonlearners,1),'linewidth',1.5,'color',colorgroup(2,:))

    ylim([0 0.05])
    xlim([0 T+1])
    xlabel('Trial','fontname','time')
    ylabel('distance [m]','fontname','time')
    % title(titles{i})

    legend(s,{'Learners','Non Learners'})
    legend boxoff
    set(gca, 'FontName', 'times');

    set(gcf,'pos',[300 0 150 190])

    %%%%%%%%%%%%%%%%%%%%% Action Analysis
    % p(a=1) collab; p(a=2) collab; p(a=0)
    purple = [138,43,226]./255;
    orange = [230, 149, 0]./255;
    light_green = [0.1 0.9 0.1];
    green = [0.1 0.6 0.1];
    light_blue = [51 153 255]./255;
    blue = [0 102 204]./255;
    alpha_shade = 0.3;
    lw = 1.5;
    trials_ep = 10;
    n_ep = T/trials_ep;

    % Actions over time
    for d = 1:n_dyads
    for s = 1:n_sim
        figure
        set(gcf,'pos',[0 0 190 110])
        hold on
        stairs(ACTION_pl1{d}(s,:),'color','b','linewidth',lw)
        stairs(ACTION_pl2{d}(s,:)+3,'color','r','linewidth',lw)
        yticks([0 1 2 3 4 5])
        ylim([-0.1 5.1])
        yticklabels({'11' '12' '21' '22' '12' '21'})
        title(['Actions D' num2str(d)])
        xlabel('Trials')
        xlim([0 170])
    end
end

% Actions at dyad level
for d = 1:n_dyads
    for s = 1:n_sim
        for t = 1:T
            actEE{d}(s,t) = (ACTION_pl1{d}(s,t) == 1) && (ACTION_pl2{d}(s,t) == 1);
            actLL{d}(s,t) = (ACTION_pl1{d}(s,t) == 2) && (ACTION_pl2{d}(s,t) == 2);

            actEL{d}(s,t) = (ACTION_pl1{d}(s,t) == 1) && (ACTION_pl2{d}(s,t) == 2);
            actLE{d}(s,t) = (ACTION_pl1{d}(s,t) == 2) && (ACTION_pl2{d}(s,t) == 1);

            actEM{d}(s,t) = (ACTION_pl1{d}(s,t) == 1) && (ACTION_pl2{d}(s,t) == 0);
            actLM{d}(s,t) = (ACTION_pl1{d}(s,t) == 2) && (ACTION_pl2{d}(s,t) == 0);

            actME{d}(s,t) = (ACTION_pl1{d}(s,t) == 0) && (ACTION_pl2{d}(s,t) == 1);
            actML{d}(s,t) = (ACTION_pl1{d}(s,t) == 0) && (ACTION_pl2{d}(s,t) == 2);

            actMM{d}(s,t) = (ACTION_pl1{d}(s,t) == 0) && (ACTION_pl2{d}(s,t) == 0);

        end
        act_EE_or_LL{d}(s,:) = actEE{d}(s,:) + actLL{d}(s,:); % collab
        act_EL_or_LE{d}(s,:) = actEL{d}(s,:) + actLE{d}(s,:); % oscill
        act_EM_LM_ME_LE{d}(s,:) = actEM{d}(s,:) + actLM{d}(s,:) + actME{d}(s,:) + actML{d}(s,:) + actMM{d}(s,:); % 1collab
        act_nocollab{d}(s,:) = act_EL_or_LE{d}(s,:) + act_EM_LM_ME_LE{d}(s,:);

    end
end

% Probabilities
ntrials_ep = 10;
for d = 1:n_dyads
    for s = 1:n_sim
        for e = 1:n_ep
            currtrials = (1:ntrials_ep) + (e-1)*ntrials_ep ;
            pEEorLL{d}(s,e) = mean(act_EE_or_LL{d}(s,currtrials));
            pELorLE{d}(s,e) = mean(act_EL_or_LE{d}(s,currtrials));
            pEMorLMorMEorLE{d}(s,e) = mean(act_EM_LM_ME_LE{d}(s,currtrials));
            pMM{d}(s,e) = mean(actMM{d}(s,currtrials));
            pnocollab{d}(s,e) = mean(act_nocollab{d}(s,currtrials));

        end

        % figure
        % set(gcf,'pos',[0 0 150 150])
        % hold on
        % stairs(pEEorLL{d}(s,:),'Color',orange,'LineWidth',lw)
        % stairs(pELorLE{d}(s,:),'Color','magenta','LineWidth',lw)
        % stairs(pEMorLMorMEorLE{d}(s,:),'Color','black','LineWidth',lw)
        % %stairs(pMM(d,:),'Color','black','LineWidth',lw)
        % title(['P(diadic actions) D' num2str(s)])
        % xlabel('Epochs')
        % ylabel('Probability')
        % xlim([0 17])
        % ylim([0 1])
        % legend('collab','cycling','not collab')
        % legend boxoff
    end
end

% % Population figure
% for d = 1:n_groups
%     figure
%     set(gcf,'pos',[0 0 150 150])
%     hold on
%     patch([1:n_ep n_ep:-1:1],[mean(pEEorLL{d})-std(pEEorLL{d})./sqrt(n_sim) fliplr(mean(pEEorLL{d})+std(pEEorLL{d})./sqrt(n_sim))],...
%         orange,'FaceAlpha',alpha_shade,'EdgeAlpha',0)
%     patch([1:n_ep n_ep:-1:1],[mean(pELorLE{d})-std(pELorLE{d})./sqrt(n_sim) fliplr(mean(pELorLE{d}) + std(pELorLE{d})./sqrt(n_sim))],...
%         'magenta','Facealpha',alpha_shade,'EdgeAlpha',0)
%     patch([1:n_ep n_ep:-1:1],[mean(pEMorLMorMEorLE{d}) - std(pEMorLMorMEorLE{d})./sqrt(n_sim) fliplr(mean(pEMorLMorMEorLE{d}) + std(pEMorLMorMEorLE{d})./sqrt(n_sim))],...
%         'k','Facealpha',alpha_shade,'EdgeAlpha',0)
% 
%     p(1) = plot(mean(pEEorLL{d}),'Color',orange,'LineWidth',lw)
%     p(2) = plot(mean(pELorLE{d}),'Color','magenta','LineWidth',lw)
%     p(3) = plot(mean(pEMorLMorMEorLE{d}),'Color','black','LineWidth',lw)
% 
%     title(['P(diadic actions)'])
%     xlabel('Epochs')
%     xlim([0 17])
%     ylim([0 1])
%     legend(p,{'collab','cycling','not collab'})
%     legend boxoff
% end

% Population figure -- Learners DYADS!!
learnersdyads = intersect(learners1,learners2);
ndyads_learners = size(learnersdyads,1);
ndyads_nonlearners = 18-ndyads_learners;
n_sub = [ndyads_learners ndyads_nonlearners];

nonlearnersdyads = 1:18;
nonlearnersdyads(learnersdyads) = [];

pEEorLL_mat = reshape(cell2mat(pEEorLL),17,18)';
pELorLE_mat = reshape(cell2mat(pELorLE),17,18)';
pnocollab_mat = reshape(cell2mat(pnocollab),17,18)';

pEEorLL_new{1} = pEEorLL_mat(learnersdyads,:);
pEEorLL_new{2} = pEEorLL_mat(nonlearnersdyads,:);
pnocollab_new{1} = pnocollab_mat(learnersdyads,:);
pnocollab_new{2} = pnocollab_mat(nonlearnersdyads,:);

tits = {'Learners','Non Learners'}
LEARNcol =[255 102 102]./255; %[0 1 0];
NOLEARNcol = [59 179 0]./255;
colcol = [LEARNcol; NOLEARNcol];

figure
for g = 1:n_groups
    
    set(gcf,'pos',[0 0 150 150])
    hold on
    patch([1:n_ep n_ep:-1:1],[mean(pEEorLL_new{g})-std(pEEorLL_new{g})./sqrt(n_sub(g)) fliplr(mean(pEEorLL_new{g})+std(pEEorLL_new{g})./sqrt(n_sub(g)))],...
         colcol(g,:),'FaceAlpha',alpha_shade,'EdgeAlpha',0)

    % % patch([1:n_ep n_ep:-1:1],[mean(pnocollab_new{g}) - std(pnocollab_new{g})./sqrt(n_sub(g)) fliplr(mean(pnocollab_new{g}) + std(pnocollab_new{g})./sqrt(n_sub(g)))],...
    % %     light_green,'Facealpha',alpha_shade,'EdgeAlpha',0)

    p(1) = plot(mean(pEEorLL_new{g}),'Color',colcol(g,:),'LineWidth',lw)
    % p(2) = plot(mean(pnocollab_new{g}),'Color',light_green,'LineWidth',lw)

    % title(tits{g})
    xlabel('Epochs')
    ylabel('Probability')
    xlim([1 17])
    ylim([0 1])
    % legend(p,{'Collaborative','Solo'})
    legend boxoff
    set(gca,'fontname','times')

end

% Stat
mdlearners_final = mean(mdlearners(:,(end-9):end),2);
mdnonlearners_final = mean(mdnonlearners(:,(end-9):end),2);



end







