%function [c1,c2] = get_corr2(nn,u1,u1hat,u2,u2hat)
function [c1,c2,c1x,c1y,c2x,c2y] = get_corr2(nn,u1,u1hat,u2,u2hat)

T = size(u1,2);
n_nodes = length(u1(:,1))/2;

x_nodes = 1:2:(n_nodes*2);
y_nodes = 2:2:(n_nodes*2);

for t=1:T
    c1temp = corrcoef(u1(:,t),u1hat(:,t));
    c1(t) = c1temp(1,2).^2;
    c2temp = corrcoef(u2(:,t),u2hat(:,t));
    c2(t) = c2temp(1,2).^2;

    c1x_tmp = corrcoef(u1(x_nodes,t),u1hat(x_nodes,t));
    c1x(t) = c1x_tmp(1,2).^2;
    c1y_tmp = corrcoef(u1(y_nodes,t),u1hat(y_nodes,t));
    c1y(t) = c1y_tmp(1,2).^2;

    c2x_tmp = corrcoef(u2(x_nodes,t),u2hat(x_nodes,t));
    c2x(t) = c2x_tmp(1,2).^2;
    c2y_tmp = corrcoef(u2(y_nodes,t),u2hat(y_nodes,t));
    c2y(t) = c2y_tmp(1,2).^2;

end

