function t = task(Rii,Ri_i,R_i_i,ri,r_i,z)
% TASK Creates a static quadratic cost function
%  
% Ji(u1,u2) = 1/2 [ ui'*Rii*ui + u-i'*Ri-i*ui + u-i'*R-i-i*u-i ] +
% ri*ui+r-i*u-i + zi
%

nt = length(z);


if nt==1 % just a quadratic game
    t.Rii = Rii;
    t.Ri_i = Ri_i;
    t.R_i_i = R_i_i;
    t.ri = ri;
    t.r_i = r_i;
    t.z = z;
    t.usize = length(Rii);
    t.xsize = length(R_i_i);

    detRii = det(Rii);

    if detRii>0
        t.type = 1; %the cost function is for coordination 
    else
        t.type = 2; %the cost function is for competition 
    end
    
else
   for n=1:nt
      t.Rii{n}=Rii{n};
      t.Ri_i{n} = Ri_i{n};
      t.R_i_i{n} = R_i_i{n};
      t.ri{n} = ri{n};
      t.r_i{n} = r_i{n};
      t.z{n} = z{n}; 
      detRii{n} = det(Rii{n});
   end
   t.usize = length(Rii{1});
   t.xsize = length(R_i_i{1});
   if detRii{1}>0
        t.type = 1;
    else
        t.type = 2;
   end
   
end

t.nt = nt;
t = class(t,'task');

end

