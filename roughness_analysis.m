Fig = openfig('vanderpol.fig');
axis equal
hold on 
grid on

% Plot the synthesized optimal controlled trajectory
plot_trajectory_optimal = scatter(trajectory_optimal(1,:),trajectory_optimal(2,:),10, 'r','filled');
% plot_trajectory_optimal = plot(trajectory_optimal(1,:),trajectory_optimal(2,:),'r','LineWidth',1);
% pause()

%% Apply discretization
dt = 0.001;
x0 = trajectory_optimal(:,1);
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

% Plot the discretized trajectory 
plot_discretized_trajectory = scatter(x_traj(1,:),x_traj(2,:),6, 'g','filled');


