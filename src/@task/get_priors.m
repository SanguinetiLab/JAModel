function prior = get_priors(task,x,P,lambda)

Mu = get_center(task,x);

if task.nt==1
    prior = 1;
else
    for n=1:task.nt
        q(n) = -0.5*Mu(:,n)'*task.Rii{n}*Mu(:,n) + ...
            task.r_i{n}*x + task.z{n} + 0.5*trace(task.R_i_i{n}.*P) + ...
            0.5*x'*task.R_i_i{n}*x; %ok
        zz(n)=-0.5*log(det(task.Rii{n}))-q(n)./lambda;
    end
    switch task.nt
    case 2
        prior(1) = 1/(1+exp(zz(2)-zz(1)));
        prior(2) = 1/(1+exp(zz(1)-zz(2)));
    otherwise
        %pr = exp(zz);
        %pr = pr*10;
        %prior = pr/sum(pr);
        for n=1:task.nt
            prior(n) = 1./sum(exp(zz-zz(n))); %ok Pr(si = SKi|xi)
        end
end
end

prior(prior<10*eps)=10*eps;

    