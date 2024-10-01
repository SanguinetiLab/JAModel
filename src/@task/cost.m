function Ji = cost(t,ui,u_i)
% COST calculates the cost of a task (part of a quadratic game) for given actions U1,U2
if t.nt==1
    Ji(1) = 1/2*(ui'*t.Rii*ui + 2*u_i'*t.Ri_i*ui + u_i'*t.R_i_i*u_i) + ...
        t.ri*ui + t.r_i*u_i + t.z;
else
    for n=1:t.nt
        Ji{n} = 1/2*(ui'*t.Rii{n}*ui + 2*u_i'*t.Ri_i{n}*ui + u_i'*t.R_i_i{n}*u_i)...
            + t.ri{n}*ui + t.r_i{n}*u_i + t.z{n};

        R=1/2*(ui'*t.Rii{n}*ui + 2*u_i'*t.Ri_i{n}*ui + u_i'*t.R_i_i{n}*u_i)
        ri=t.ri{n}*ui
        r_i=t.r_i{n}*u_i
        zi= t.z{n}
        R+ri+r_i+zi
    end
end

end

