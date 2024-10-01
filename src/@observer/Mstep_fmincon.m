function e = Mstep_fmincon(e,c,t,K,u,u_other)

p0 = ;
p_low = ;
p_up = ;

[p_new,fval] = fmincon(@ll_decisiomaking,p0,[],[],[],[],p_low,p_up,...
    [],[],e,c,t,K,u,u_other);


end

