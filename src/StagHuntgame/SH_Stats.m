%% Stat comparison of the different simulations
clear all
clc
close all

n_sims = 4
n_trials = 40;
n_blocks = 50;
n_dyads = 8;

load sim3.mat
for d = 1:n_dyads
    for b = 1:n_blocks
        b_trials =((b-1)*n_trials+1):b*n_trials;
        n_HH_block(b) = length(find(HH(d,b_trials)));
        n_SS_block(b) = length(find(SS(d,b_trials)));
    end
    n_HH(1,d) = mean(n_HH_block);
    n_SS(1,d) = mean(n_SS_block);
end

clear b_trials n_HH_block n_SS_block
load sim5.mat
for d = 1:n_dyads
    for b = 1:n_blocks
        b_trials =((b-1)*n_trials+1):b*n_trials;
        n_HH_block(b) = length(find(HH(d,b_trials)));
        n_SS_block(b) = length(find(SS(d,b_trials)));
    end
    n_HH(2,d) = mean(n_HH_block);
    n_SS(2,d) = mean(n_SS_block);
end

clear b_trials n_HH_block n_SS_block
load sim7.mat
for d = 1:n_dyads
    for b = 1:n_blocks
        b_trials =((b-1)*n_trials+1):b*n_trials;
        n_HH_block(b) = length(find(HH(d,b_trials)));
        n_SS_block(b) = length(find(SS(d,b_trials)));
    end
    n_HH(3,d) = mean(n_HH_block);
    n_SS(3,d) = mean(n_SS_block);
end

clear b_trials n_HH_block n_SS_block
load sim9.mat
for d = 1:n_dyads
    for b = 1:n_blocks
        b_trials =((b-1)*n_trials+1):b*n_trials;
        n_HH_block(b) = length(find(HH(d,b_trials)));
        n_SS_block(b) = length(find(SS(d,b_trials)));
    end
    n_HH(4,d) = mean(n_HH_block);
    n_SS(4,d) = mean(n_SS_block);
end

p = kruskalwallis(n_HH');
p = kruskalwallis(n_SS');

p = ranksum(n_SS(1,:)',n_HH(1,:)');
p = ranksum(n_SS(2,:)',n_HH(2,:)');
p = ranksum(n_SS(3,:)',n_HH(3,:)');
p = ranksum(n_SS(4,:)',n_HH(4,:)');

%% SIGMAY EFFECT
p = ranksum([n_SS(1,:)'; n_SS(2,:)'],[n_SS(3,:)'; n_SS(4,:)']);
[h,p] = kstest(n_SS(1,:)); [h,p] = kstest(n_SS(2,:)); [h,p] = kstest(n_SS(3,:)); [h,p] = kstest(n_SS(4,:));
anova2([[n_SS(1,:)'; n_SS(2,:)'],[n_SS(3,:)'; n_SS(4,:)']],8)

%% SIGMAx EFFECT
p = ranksum([n_SS(1,:)'; n_SS(3,:)'],[n_SS(2,:)'; n_SS(4,:)']);
anova2([[n_SS(1,:)'; n_SS(3,:)'],[n_SS(2,:)'; n_SS(4,:)']],8)










