function [H,L,sigmay,A,sigmax,P0,x0,lambda,decay_rate] = define_player_model_from_fitting(Params,xsz,duration,nodes_num,final,start,VP)

Params.SigmaY
Params.SigmaX
Params.P0

if isfield(Params,'Astd')
    % SENSORY SYSTEM
    % sensory input equals partner action + own action + noise
    H = eye(xsz);
    L = -eye(xsz);
    sigmay = param_init_gaussian(Params.SigmaY,Params.SigmaYstd,[0 inf])*eye(xsz);%mean(SigmaY1(:,dyads))*eye(xsz1);

    % PARTNER MODEL
    % Retention rate
    A = param_init_gaussian(Params.A,Params.Astd,[0 inf])*eye(xsz);%mean(Mem1(:,dyads))*eye(xsz1);

    % Process noise
    sigmax =param_init_gaussian(Params.SigmaX,Params.SigmaXstd,[0 inf])*eye(xsz);%mean(SigmaX1(:,dyads))*eye(xsz1);

    % Initial state mean and variance
    P0 = param_init_gaussian(Params.P0,Params.P0std,[0 inf])*eye(xsz);

    dt = duration/(nodes_num+2);
    tt = 0:dt:duration;
    tt_half = 0:dt:(duration/2);

    x0 = eps.*ones(xsz,1);

    % H: the other is doing a straight line from start to end point
    for n = 1:nodes_num
        xmin_jeark(n) =  start(1) + (final(1)-start(1))*minimum_jerk_traj(tt(n+1)/duration);
        x0(n*2-1) =  xmin_jeark(n) ;
    end
    % H: the other is doing my task
    for n = 1:(nodes_num/2)
        ymin_jeark1(n) = start(2) + (VP(2)-start(2))*minimum_jerk_traj(tt_half(n+1)/(duration/2));
        x0([n*2 (2*nodes_num-2*(n-1))]) = ymin_jeark1(n);
    end

    % ACTION CONTROLLER
    % Temperature d0ecreases during experiment
    lambda = param_init_gaussian(Params.lambda0,Params.lambda0std,[0 inf]);

    decay_rate = param_init_gaussian(Params.a,Params.astd,[0 inf]);
else
    % SENSORY SYSTEM
    % sensory input equals partner action + own action + noise
    H = eye(xsz);
    L = -eye(xsz);
    sigmay = ((Params.SigmaY)*eye(xsz));%mean(SigmaY1(:,dyads))*eye(xsz1);

    % PARTNER MODEL
    % Retention rate
    A = Params.A*eye(xsz);%mean(Mem1(:,dyads))*eye(xsz1);

    % Process noise
    sigmax =((Params.SigmaX)*eye(xsz));%mean(SigmaX1(:,dyads))*eye(xsz1);

    % Initial state mean and variance
    P0 = ((Params.P0)*eye(xsz));

    dt = duration/(nodes_num+2);
    tt = 0:dt:duration;
    tt_half = 0:dt:(duration/2);

    x0 = eps.*ones(xsz,1);

    % H: the other is doing a straight line from start to end point
    for n = 1:nodes_num
        xmin_jeark(n) =  start(1) + (final(1)-start(1))*minimum_jerk_traj(tt(n+1)/duration);
        x0(n*2-1) =  xmin_jeark(n) ;
    end
    % H: the other is doing my task
    for n = 1:(nodes_num/2)
        ymin_jeark1(n) = start(2) + (VP(2)-start(2))*minimum_jerk_traj(tt_half(n+1)/(duration/2));
        x0([n*2 (2*nodes_num-2*(n-1))]) = ymin_jeark1(n);
    end

    % ACTION CONTROLLER
    % Temperature d0ecreases during experiment
    lambda = Params.lambda0;

    decay_rate = Params.a;

end
