function [posx,posy] = get_position_np(traj,interval,hand,n)
% Get positions normalized based of path length
% Cecilia 2019
% n = number of point for the undersampling

% Path Normalization
temp_path = 0;
path(1) = 0;
if traj.trsize == 2
    for i = 2:length(interval)
        temp_path = sqrt((diff(traj.pos(interval(i-1:i),2*(hand-1)+1))).^2+(diff(traj.pos(interval(i-1:i),2*(hand-1)+2))).^2);
        path(i) = path(i-1) + temp_path;
    end
else  
    for i = 2:length(interval)
        temp_path = sqrt((diff(traj.pos(interval(i-1:i),3*(hand-1)+2))).^2+(diff(traj.pos(interval(i-1:i),3*(hand-1)+3))).^2);
        path(i) = path(i-1)+temp_path;
    end
end
path = path./max(path);

% [md1,i1]= min(sqrt((traj.pos(interval,3*(hand-1)+2)-traj.viapoints(1,2)).^2+...
%                    (traj.pos(interval,3*(hand-1)+3)-traj.viapoints(1,3)).^2));
% [md2,i2]= min(sqrt((traj.pos(interval,3*(hand-1)+2)-traj.viapoints(2,2)).^2+...
%                    (traj.pos(interval,3*(hand-1)+3)-traj.viapoints(2,3)).^2));             
% % tc1 = traj.time(interval(i1))-traj.time(interval(1));
% % tc2 = traj.time(interval(i2))-traj.time(interval(1));
% time = traj.time(interval)-traj.time(interval(1));
% time = time./max(time);
% tc1 = time(i1);
% tc2 = time(i2);
% half_path_ind = find(path>=0.5,1);

% figure
% set(gcf,'pos',[0 0 400 1000])
% subplot(3,1,1)
% plot(time,path)
% hold on
% plot(time(half_path_ind),path(half_path_ind),'*')
% xline(tc1,'-',{'tc1'})
% xline(tc2,'-',{'tc2'})
% box off
% xlabel('Normalized Time')
% ylabel('Normalized Path')
% title('Path/Time')
% 
% subplot(3,1,2)
% plot(traj.pos(interval,3*(hand-1)+2),traj.pos(interval,3*(hand-1)+3))
% % plot(time,traj.pos(interval,3*(hand-1)+3))
% box off
% title('Trajectory')
% speed = sqrt(traj.vel(interval,3*(hand-1)+2).^2+traj.vel(interval,3*(hand-1)+3).^2);
% 
% subplot(313)
% plot(time,speed)
% xlabel('Normalized Time')
% ylabel('Speed (m/s)')
% title('Speed')

% undersamples = 0:0.01:1; %% ISSUE: LA TRAIETTORIA MEDIA NON PASSA PER I
% VP
undersamples = linspace(0,1,n);

for i =1:1:length(undersamples)
    idx(i) = find(path>=undersamples(i),1);
    if i == length(undersamples)
        idx(i) = length(path);
    end
end

% %% tentativo allineare le traiettorie in modo da avere traiettoria media che passa per i via points
% undersamples = linspace(0,0.1,100);
% for i =1:length(undersamples)
%     temp = find(traj.pos(interval,3*(hand-1)+2)>=undersamples(i),1);
%     if size(temp)>1
%         diff_temp = temp - idx(i-1);
%         [m,ind_m] = min(diff_temp);
%         idx(i) = temp(ind_m);
%     elseif size(temp) == 1
%         idx(i) = temp;
%     end
%     if i == length(undersamples)
%         idx(i) = length(traj.pos(interval,3*(hand-1)+2));
%     end
% end
% for i = 1:length(idx)
%     if idx(i) == 0
%         idx(i) = idx(i-1);
%     end
% end
% 


switch traj.trsize
    case 3, %3D
        tr = traj.pos(interval,3*(hand-1)+2:3*(hand-1)+3);
        posx = tr(idx,1); 
        posy = tr(idx,2);
%         under_tr = [posx posy];
    case 2,  %Planar
        tr = traj.pos(interval,1:2); %% Da vedere meglio
        posx = tr(idx,1);
        posy = tr(idx,2);
%         under_tr = [posx posy];
end 



