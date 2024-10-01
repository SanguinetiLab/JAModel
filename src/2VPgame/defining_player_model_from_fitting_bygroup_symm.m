function Params = defining_player_model_from_fitting_bygroup_symm(datadir,suffix,n_dyads,n_groups,dyads)

for d = 1:n_dyads
    data_param = load(fullfile(datadir,['Dyads_' num2str(d) '_' suffix '.mat']));
    A1(d) = data_param.Pl1_Params.A; A2(d) = data_param.Pl2_Params.A;
    SigmaY1(d) = data_param.Pl1_Params.SigmaY; SigmaY2(d) = data_param.Pl2_Params.SigmaY;
    SigmaX1(d) = data_param.Pl1_Params.SigmaX; SigmaX2(d) = data_param.Pl2_Params.SigmaX;
    P01(d) = data_param.Pl1_Params.P0; P02(d) = data_param.Pl2_Params.P0;
    lambda01(d) = data_param.Pl1_Params.lambda0; lambda02(d) = data_param.Pl2_Params.lambda0;
    a1(d) = data_param.Pl1_Params.a; a2(d) = data_param.Pl2_Params.a;
end

for g = 1:n_groups
    
    curr_dyads = dyads{g};
    for d = curr_dyads
        pp.A = mean([A1(curr_dyads) A2(curr_dyads)]);
        pp.Astd = std([A1(curr_dyads) A2(curr_dyads)]);

        pp.SigmaY = mean([SigmaY1(curr_dyads) SigmaY2(curr_dyads)]);
        pp.SigmaYstd = std([SigmaY1(curr_dyads) SigmaY2(curr_dyads)]);

        pp.SigmaX = mean([SigmaX1(curr_dyads) SigmaX2(curr_dyads)]);
        pp.SigmaXstd = std([SigmaX1(curr_dyads) SigmaX2(curr_dyads)]);

        pp.P0 = mean([P01(curr_dyads) P02(curr_dyads)]);
        pp.P0std = std([P01(curr_dyads) P02(curr_dyads)]);

        pp.lambda0 = mean([lambda01(curr_dyads) lambda02(curr_dyads)]);
        pp.lambda0std = std([lambda01(curr_dyads) lambda02(curr_dyads)]);

        pp.a = mean([a1(curr_dyads) a2(curr_dyads)]);
        pp.astd = std([a1(curr_dyads) a2(curr_dyads)]);

    end

    for i = 1:length(curr_dyads)
        dd = curr_dyads(i);
        Params{dd} = pp;
    end

end



