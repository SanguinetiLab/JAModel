function [converges] = convergence(game)
% CONVERGENCE(GAME) Determines convergence of a static quadratic game through FP
% For each NE of the game, determines if the FP converges to that NE


% Ji(u1,u2) = 1/2 [ ui'*Rii*ui + 2 u-i'*Ri-i*ui + u-i'*R-i-i*u-i ] +
% ri*ui+r-i*u-i + zi
%
% reaction line: dJi/dui = Rii*ui + R-i*u-i +ri = 0

% Method: for a quadratic game, FP can be seen as a Jacobi iterative method

% Get the game parameters
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

if nt1 ~= nt2
    error('inconsistent number of NEs!')
else
    if nt1==1 % one NE
        R = [Rii_1 Ri_i_1; Ri_i_2 Rii_2];

        % Decompose R in a (block) diagonal and a non-diagonal part:
        D = [Rii_1 0*Ri_i_1; 0*Ri_i_2 Rii_2];
        ND = [0*Rii_1 Ri_i_1; Ri_i_2 0*Rii_2];
        C = inv(D)*ND;
        eigvals = eig(C);
        converges = all(norm(eigvals))<1

        r = [ri_1';ri_2'];
        u_ne = -inv(R)*r;
        u1ne = u_ne(1:u1sz);
        u2ne = u_ne((u1sz+1):end);

    else
        for n1=1:nt1
            R = [Rii_1{n1} Ri_i_1{n1}; Ri_i_2{n1} Rii_2{n1}];
            
            % Decompose R in a (block) diagonal and a non-diagonal part:
            D = [Rii_1{n1} 0*Ri_i_1{n1}; 0*Ri_i_2{n1} Rii_2{n1}];
            ND = [0*Rii_1{n1} Ri_i_1{n1}; Ri_i_2{n1} 0*Rii_2{n1}];
            C = inv(D)*ND;
            eigvals = eig(C);
            converges(n1) = all(abs(eigvals)<1)
   
            r = [ri_1{n1}';ri_2{n1}'];
            u_ne = -inv(R)*r;
            u1ne{n1,n1} = u_ne(1:u1sz);
            u2ne{n1,n1} = u_ne((u1sz+1):end);
        end
    end
end

