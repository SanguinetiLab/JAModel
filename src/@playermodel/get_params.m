function pars = get_params(model)
% get model parameters
for p=1:length(model.parnames)
   ['pars{p}= model.',parnames{p},';']
   eval(['pars{p}= model.',parnames{p},';']);
end
