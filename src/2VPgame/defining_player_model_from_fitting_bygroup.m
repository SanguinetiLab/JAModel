function Params = defining_player_model_from_fitting_bygroup(datadir,suffix,n_dyads,n_groups,n_dyads_group,VP)

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
    
    curr_dyads = (g-1)*n_dyads_group+(1:n_dyads_group);
    for d = curr_dyads
        pp.A = mean([A1(curr_dyads) A2(curr_dyads)]);
        pp.SigmaY = mean([SigmaY1(curr_dyads) SigmaY2(curr_dyads)]);
        pp.SigmaX = mean([SigmaX1(curr_dyads) SigmaX2(curr_dyads)]);
        pp.P0 = mean([P01(curr_dyads) P02(curr_dyads)]);
        pp.lambda0 = mean([lambda01(curr_dyads) lambda02(curr_dyads)]);
        pp.a = mean([a1(curr_dyads) a2(curr_dyads)]);
        Params{d} = pp;
    end
    

end



