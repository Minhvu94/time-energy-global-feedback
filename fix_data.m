clear
clc

load('u_optimal.mat','store_u_optimal')
load('T_optimal.mat','store_T_optimal')
load('trajectory_optimal.mat','store_trajectory_optimal')
%%
store_u_optimal(:,33*8+1:33*25) = mid_store_u_optimal;
store_T_optimal(:,33*8+1:33*25) = mid_store_T_optimal;
store_trajectory_optimal(1,33*8+1:33*25) = mid_store_trajectory_optimal;


%% Remove the optimal data corresponding to a particular index
index = boolean([zeros(1,544),1,zeros(1,544)]);
store_trajectory_optimal(index) = [];
store_T_optimal(:,545) = [];
store_u_optimal(:,545) = [];

%% Fixing the optimal data corresponding to a particular x0/index
store_trajectory_optimal{421} = [x0,x_traj];
store_T_optimal(:,421) = T;
store_u_optimal(:,421) = u;

%% Save
save('u_optimal.mat','store_u_optimal')
save('T_optimal.mat','store_T_optimal')
save('trajectory_optimal.mat','store_trajectory_optimal')

%% Then validating if the fix was good 
Fig = openfig('vanderpol.fig');
axis equal
hold on 
grid on
% The optimal trajectory from data 
for i=1:49
    if i~=25
%     T = store_T_optimal(:,i);
%     u = store_u_optimal(:,i);
    x_traj = store_trajectory_optimal{i};
    data_traj = plot(x_traj(1,:),x_traj(2,:),'r','LineWidth',1);
scatter(x_traj(1,:),x_traj(2,:),15, 'k','filled');
%     pause(1)
%     delete(data_traj);
    end
end

%% The optimal trajectory from applying the control u for T seconds 
x = x_traj(:,1);    
x_traj = [];
for iter = 1:300
    [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 T(iter)],[x;u(iter)]); 
    x = x_trajJ_fine(end,:)'; 
    x = x(1:2);
    x_traj = [x_traj x];         
end
plot(x_traj(1,:),x_traj(2,:),'k','LineWidth',2);
