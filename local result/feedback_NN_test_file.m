%% Testing the trained NN 
x0 = [2.5;-3.5];
x = x0;    
x_traj = x0;
dt = 0.001;
step = 0;
while true 
u = myNeuralNetworkFunction(x);
[~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 dt],[x;u]); 
x = x_trajJ_fine(end,:)'; 
x = x(1:2);
x_traj = [x_traj x]; 
if norm(x-[0;0]) < 0.01
    break
end
step = step + 1
end
scatter(x_traj(1,:),x_traj(2,:),6, 'k','filled');