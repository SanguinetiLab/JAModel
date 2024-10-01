function c = get_center(task,x)
% GET_CENTER Determines the minimum of the cost function

if task.nt==1
    k = task.Ri_i *x + task.ri';
    c = -inv(task.Rii)*k;    
else
    for n=1:task.nt
      k = task.Ri_i{n} *x + task.ri{n}';
      c(:,n) = -inv(task.Rii{n})*k;
    end
end

