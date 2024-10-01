function model = playermodel(task, pars, parnames)
% Playermodel constructor
% A playemodel object is defined by task, observer and controller

model.task = task;
model.parnames = parnames;

for p=1:length(pars)
   %['model.',parnames{p},'= pars{p};']
   eval(['model.',parnames{p},'= pars{p};']);
end

model.observer = observer(model.A,model.B,model.C,model.D,model.SigmaX,model.SigmaY,model.x0,model.P0);
model.controller = controller(model.task,model.lambda0,model.a);

model = class(model,'playermodel');