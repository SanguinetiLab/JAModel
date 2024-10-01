function cost=costfunction(X,task,U,Y)

xsz = get(task,'xsize');

% pars = {exp(-1/X(1))*eye(xsz),[],eye(xsz),eye(xsz),X(2)*eye(xsz),X(3)*eye(xsz),X(3+(1:xsz)),X(3+xsz+1)*eye(xsz),X(3+xsz+2),exp(-1/X(3+xsz+3))};
pars = {X(1)*eye(xsz),[],eye(xsz),eye(xsz),X(2)*eye(xsz),X(3)*eye(xsz),X(3+(1:xsz)),X(3+xsz+1)*eye(xsz),X(3+xsz+2),X(3+xsz+3)};

parnames = {'A','B','C','D',...
    'SigmaX','SigmaY',...
    'x0','P0','lambda0','a'};

model = playermodel(task,pars,parnames);

[cost,~,~,~] = pseudolik(model,U',Y');