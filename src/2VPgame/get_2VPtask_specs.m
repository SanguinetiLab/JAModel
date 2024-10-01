function  [tcE,tcL,duration, yrange, start, final, VP1, VP2, tgterror,VPerror] = get_2VPtask_specs(mode)

switch mode
    case 'Asymm'
        % Relative Crossing time can be constant for all the dyads or can be adapted for
        % the single dyads. Here we hypothesize that are equal for all the
        % dyads
        tcE = 0.32;%0.35
        tcL = 0.68;%0.65

        duration = 2;

        % Workspace
        yrange = [-0.1 0.1];

        start = [-0.05 0];
        final = [0.05 0];

        VP1 = [-0.02  -0.03];
        VP2 = [0.02  0.03];

        tgterror = 0.005;
        VPerror =  0.0025;
    case 'Symm'

        tcE = 0.38;%0.35
        tcL = 0.62;%0.65

        duration = 2;

        % Workspace
        yrange = [-0.1 0.1];

        start = [-0.05 0];%[0 0];%
        final = [0.05 0];%[0.1 0];%

        VP1 = [0  -0.03];%[0.05  -0.03];%
        VP2 = [0  0.03];%[0.05  0.03];%

        tgterror = 0.005;
        VPerror =  0.0025;%0.0025;


end