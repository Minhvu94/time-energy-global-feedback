clear
clc
load('u_optimal.mat','store_u_optimal')
load('T_optimal.mat','store_T_optimal')
load('trajectory_optimal.mat','store_trajectory_optimal')

%%
Fig = openfig('vanderpol.fig');
axis equal
hold on 
grid on
for i=38:50
    i
    x_traj = store_trajectory_optimal{i};
    step2_current_traj = plot(x_traj(1,:),x_traj(2,:),'r','LineWidth',2);
    pause()
    delete(step2_current_traj);
end

%%
index = boolean([zeros(1,40),1,zeros(1,40)]);
store_trajectory_optimal(index) = [];
store_T_optimal(:,41) = [];
store_u_optimal(:,41) = [];
%%
store_trajectory_optimal{46} = [x0,x_traj];
store_T_optimal(:,46) = T;
store_u_optimal(:,46) = u;
%%
save('u_optimal.mat','store_u_optimal')
save('T_optimal.mat','store_T_optimal')
save('trajectory_optimal.mat','store_trajectory_optimal')

