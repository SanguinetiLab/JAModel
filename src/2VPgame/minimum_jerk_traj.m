





function xi=minimum_jerk_traj(tau)

if tau <= 0 
    xi = 0;
elseif tau < 1 & tau > 0
    xi    = tau.^3 .* (6.*tau.^2 - 15.*tau + 10);
else 
    xi = 1;
end
