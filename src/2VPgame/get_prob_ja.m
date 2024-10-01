function [collab,collab12,collab21,oscil,nocollab,oscil12_21,oscil21_12,bothignore,oneignore] = get_prob_ja(ACTION_pl1,ACTION_pl2)



T = length(ACTION_pl1);


collab = zeros(1,T); 
collab12= zeros(1,T); 
collab21= zeros(1,T); 
oscil= zeros(1,T); 
nocollab= zeros(1,T); 
oscil12_21= zeros(1,T); 
oscil21_12= zeros(1,T); 
bothignore= zeros(1,T); 
oneignore= zeros(1,T); 
for t = 1:T
    if ~(ACTION_pl1(t) == 0) && ~(ACTION_pl2(t) == 0) && (ACTION_pl1(t) == ACTION_pl2(t))
        collab(t) = 1;
        if (ACTION_pl1(t)==1)
            collab12(t) = 1;
        elseif (ACTION_pl1(t)==2)
            collab21(t) = 1;
        end

    elseif ~(ACTION_pl1(t) == 0) && ~(ACTION_pl2(t) == 0) && ~(ACTION_pl1(t) == ACTION_pl2(t))
        oscil(t) = 1;
        nocollab(t) = 1;
        if (ACTION_pl1(t) == 1) && (ACTION_pl2(t) == 2)
            oscil12_21(t) = 1;
        elseif (ACTION_pl1(t) == 2) && (ACTION_pl2(t) == 1)
            oscil21_12(t) = 1;
        end

    elseif (ACTION_pl1(t) == 0) && (ACTION_pl2(t) == 0)
        bothignore(t) = 1;
        nocollab(t) = 1;

    elseif ((ACTION_pl1(t) == 0) && ~(ACTION_pl2(t) == 0)) || ((ACTION_pl2(t) == 0) && ~(ACTION_pl1(t) == 0))
        oneignore(t) = 1;
        nocollab(t) = 1;
    end
end