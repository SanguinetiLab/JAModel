function [u,i,c]=generate_action(c,x,P,lambda)

nt = get(c.task,'nt');
type = get(c.task,'type');
curr_task = c.task;
% workspace limitations are actually useful also for non competitive games

if type==1
    
    Mu = get_center(c.task,x);
    Sigma = get_covariance(c.task,lambda);
    c.prior = get_priors(c.task,x,P,lambda);

    if nt==1 % one strategy
        % sample u from a MVN distribution
        u = mvnrnd(Mu,Sigma,1);
        % Cecilia's stupid way to limit workspace
        while u < c.umin || u > c.umax
            %u = mvnrnd(Mu,Sigma,1);
            u_minc = maxgibbs(c.task,x,c.umin,c.umax);
            u = rndgibbs(c.task,u_minc,x,c.umin,c.umax,lambda);
            %u = limit_gaussian_action(c.task,x,c.umin,c.umax,lambda,Mu,Sigma);
        end

        %limiting action to workspace boundaries

        i = 1;
    else % multiple strategies

        % if Sigma(:,:,1) < 1e-17
        %    Sigma(:,:,1) = 1e-17;
        %    Sigma(:,:,2) = 1e-17;
        % end
        gm = gmdistribution(Mu',Sigma,c.prior);
        % sample u from a mixture of Gaussians
        %c.prior
        %Mu


        [u,i] = random(gm);
    end
else

    u_minc = maxgibbs(c.task,x,c.umin,c.umax);
    u = rndgibbs(c.task,u_minc,x,c.umin,c.umax,lambda);
    i = 1;
    c.prior = get_priors(c.task,x,P,lambda);

end


end








