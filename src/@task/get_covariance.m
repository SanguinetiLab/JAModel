function Sigma = get_covariance(task,lambda)

if task.nt==1
    Sigma = lambda*inv(task.Rii);
    Sigma = (Sigma+Sigma')/2; % To make sure it is symmetric...
else
    for n=1:task.nt
      Sigma(:,:,n) = lambda*inv(task.Rii{n});
      Sigma(:,:,n) = (Sigma(:,:,n)+Sigma(:,:,n)')/2; % To make sure it is symmetric...
    end
end