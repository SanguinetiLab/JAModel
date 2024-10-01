function ACTION_pl = get_action_class(th_md,mdij,tci1,tci2)

T = length(mdij);

for t = 1:T
    if mdij(t) <= th_md
        if tci1(t) <= tci2(t)
            ACTION_pl(t) = 1;
        else
            ACTION_pl(t) = 2;
        end
    else
        ACTION_pl(t) = 0;
    end
end