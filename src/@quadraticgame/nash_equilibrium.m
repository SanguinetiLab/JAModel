function [u1ne,u2ne] = nash_equilibrium(game)
% NASH_EQUILIBRIUM Computes Nash equilibrium for a static quadratic game

% Ji(u1,u2) = 1/2 [ ui'*Rii*ui + 2 u-i'*Ri-i*ui + u-i'*R-i-i*u-i ] +
% ri*ui+r-i*u-i + zi
%
% dJi/dui = Rii*ui + R-i*u-i +ri = 0

Rii_1 = get(game.task1,'Rii');
Ri_i_1 = get(game.task1,'Ri_i');
ri_1 = get(game.task1,'ri');

Rii_2 = get(game.task2,'Rii');
Ri_i_2 = get(game.task2,'Ri_i');
ri_2 = get(game.task2,'ri');

u1sz = get(game.task1,'usize');
u2sz = get(game.task2,'usize');

nt1 = get(game.task1,'nt');
nt2 = get(game.task1,'nt');

if nt1==1 & nt2==1
    R = [Rii_1 Ri_i_1; Ri_i_2 Rii_2];
    r = [ri_1';ri_2'];
    u_ne = -inv(R)*r;

    u1ne = u_ne(1:u1sz);
    u2ne = u_ne((u1sz+1):end);
else
    for n1=1:nt1
        for n2 = 1:nt2
            R = [Rii_1{n1} Ri_i_1{n1}; Ri_i_2{n2} Rii_2{n2}];
            R = [Rii_1{n1} Ri_i_1{n1}; Ri_i_2{n2} Rii_2{n2}];

            r = [ri_1{n1}';ri_2{n2}'];
            u_ne = -inv(R)*r;
            u1ne{n1,n2} = u_ne(1:u1sz);
            u2ne{n1,n2} = u_ne((u1sz+1):end);
        end
    end
end

