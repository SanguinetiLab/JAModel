function g = quadraticgame(task1,task2)
% QUADRATICGAME Creates a static quadratic game

% TODO: must check that the dimensions are compatible: task1.usize must be equal
% to task2.xsize and vice versa

g.task1 = task1;
g.task2 = task2;

g = class(g,'quadraticgame');

end

