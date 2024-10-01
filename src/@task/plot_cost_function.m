function [J,react_curve] = plot_cost_function(task,action_space)

actions = action_space;%[action_range(1):action_range(2)];

nt = get(task,'nt');

[ui,u_i] = meshgrid(actions,actions);
for i = 1:size(ui,1)
    for j = 1:size(u_i,2)
        Ji_temp = cost(task,ui(i,j),u_i(i,j));
        for d = 1:nt
            if nt == 1
                eval(['Ji_ch' num2str(d) '(i,j) = Ji_temp(d);'])
            else
                eval(['Ji_ch' num2str(d) '(i,j) = Ji_temp{d};'])
            end
            eval(['J{d} = Ji_ch' num2str(d) ';']) 
        end
    end
end

    %% Cost Function
    figure
    for d = 1:nt
        subplot(1,nt,d)
        eval(['surf(ui,u_i,J{' num2str(d) '})'])
        tit = ['Ji - NE ' num2str(d)];
        xlabel('u_i')
        ylabel('u_{-i}')
        title(tit)
        axis square
        eval(['jmax(d) = max(max(Ji_ch' num2str(d) '));'])
        view(2)  
    end

    for d = 1:nt
        clim([0 max(jmax)])
    end




    %% Reaction Curves
    figure

    for d = 1:nt
        
        subplot(1,nt,d)
        for i = 1:length(actions)
            eval(['[mn,idxmin] = (min(Ji_ch' num2str(d) '(i,:)));']);
            reaction_curve(i) = actions(idxmin);
        end
        plot(actions,reaction_curve)

        react_curve{nt} = reaction_curve;
    end





