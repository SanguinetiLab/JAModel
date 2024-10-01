function  Pl_Params = struct_params(parnames_savings,nodes_num,X)

prec_idx = 0;
for p = 1:length(parnames_savings)
    idx = prec_idx + 1;

    switch parnames_savings{p}
        case 'x0'
            idx = (prec_idx + 1):(prec_idx(end) + nodes_num*2);
    end

    eval(['Pl_Params.' parnames_savings{p} ' = X(idx);'])

    prec_idx = idx(end);
end