%% Open the Neural Network Start GUI with this command:  nnstart

%%
clear
clc
% myCell1 = cell(1, 10);  
% myCell2 = {5,5,5,5,5};
% myCell1(1,1:5) = myCell2;
%%
load('u_optimal.mat','store_u_optimal')
load('T_optimal.mat','store_T_optimal')
load('trajectory_optimal.mat','store_trajectory_optimal')

%%
% Fig = openfig('vanderpol.fig');
axis equal
hold on 
grid on
for i=33*25:33*33-1
    
%     i
    x_traj = store_trajectory_optimal{i};
%     scatter(x_traj(1,:),x_traj(2,:),10, 'r','filled');
    step2_current_traj = plot(x_traj(1,:),x_traj(2,:),'r','LineWidth',1);
%     pause()
%     delete(step2_current_traj);
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

%% Task 1: Prepare training data via finer discretization 
Fig = openfig('vanderpol.fig');
axis equal
hold on 
grid on

training_X = [];
training_U = [];
for episode = 1:49
if episode~=25
% Plot the synthesized optimal controlled trajectory
trajectory_optimal = store_trajectory_optimal{episode};
% plot_trajectory_optimal = scatter(trajectory_optimal(1,:),trajectory_optimal(2,:),15, 'r','filled');
% plot_trajectory_optimal = plot(trajectory_optimal(1,:),trajectory_optimal(2,:),'r','LineWidth',1);
% pause()

% Apply discretization
u = store_u_optimal(:,episode);
T = store_T_optimal(:,episode);
dt = 0.005;

x0 = trajectory_optimal(:,1);
% Print initial state being discretized 
data = [episode, x0(1), x0(2)];
formatSpec = 'episode = %.0f; x = %.2f; y = %.2f \n';
fprintf(formatSpec, data)

x = x0;    
x_traj = x0;
u_traj = [];
for iter = 1:300
    % Divide a large time inverval into smaller ones, each approxi. dt
    if T(iter) > dt
        incremental_num = floor(T(iter)/dt);
        for ii = 1:incremental_num
            [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 T(iter)/incremental_num],[x;u(iter)]); 
            x = x_trajJ_fine(end,:)'; 
            x = x(1:2);
            x_traj = [x_traj x];  
            u_traj = [u_traj; u(iter)];
        end
    else
        [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 T(iter)],[x;u(iter)]); 
        x = x_trajJ_fine(end,:)'; 
        x = x(1:2);
        x_traj = [x_traj x]; 
        u_traj = [u_traj; u(iter)];
    end
end
training_X = [training_X x_traj(:,1:end-1)];
training_U = [training_U u_traj'];

% Plot the discretized trajectory 
plot_discretized_trajectory = scatter(x_traj(1,:),x_traj(2,:),6, 'g','filled');
% plot_trajectory_optimal = scatter(trajectory_optimal(1,:),trajectory_optimal(2,:),15, 'r','filled');
pause(0.01)
% 
% delete(plot_trajectory_optimal)
% delete(plot_discretized_trajectory)
end
end
%% Testing the trained NN 
Fig = openfig('vanderpol.fig');
grid on 

[X,Y] = meshgrid(-1:0.5:1 , -1:0.5:1);
grids = [X(:)  Y(:)];

x0 = [3;-3];
% for rep=1:length(grids)
% x0 = grids(rep,:)';
x = x0;    
x_traj = x0;
dt = 0.005;
step = 0;
for step=1:10000
u = myNeuralNetworkFunction(x);
[~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 dt],[x;u]); 
x = x_trajJ_fine(end,:)'; 
x = x(1:2);
x_traj = [x_traj x]; 
if norm(x-[0;0]) < 0.01
    break
end
end
step
scatter(x_traj(1,:),x_traj(2,:),6, 'g','filled');
pause(0.01)

% end



















