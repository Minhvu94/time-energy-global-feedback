clear
clc
close all
syms xi1 xi2 xi3 dt

F = [xi2;
     -xi1+xi2-xi1^2*xi2+xi3;
      0];
  
F = simplify(F);
     
p = 6;
[Phi,Psi_p,JPhi,dPhi_dh_syms] = compute_Phi_and_JPhi2(p,F,[xi1 xi2 xi3],dt);
%%
[X,Y] = meshgrid(1:0.5:10,0:20);
Z = X.*Y;
surf(X,Y,Z)
hold on 
plot3([10,1],[0,20],[0,20])
%% Part 1:
tic
dt = 0.01;
T = 5;
iter_max = ceil(T/dt);
x0 = [-3;3];
% Using orignal iterative:
%     (-2.00;3.00); U = 27.7694; dt = 0.01; T = 4.00; xn-xf = 0.0000;  --> U^2*dt=7.7114;
%     (-2.00;3.00); U = 13.9331; dt = 0.02; T = 8.00; xn-xf = 0.0000;  --> U^2*dt=3.8826; (approach the further manifold) 
% Adjusting dt:
% U = 27.7692; dt = 0.01; ddt = 0.000000; T = 4.000 ; xn-xf = 0.0000; improve = 0.0000001; --> U^2*dt=7.7113;
% U = 11.7078; dt = 0.02; ddt = 0.000000; T = 8.000 ; xn-xf = 0.0000; improve = 0.0000001; --> U^2*dt=2.7415; (stay at the closer manifold) 
% U = 8.2058 ; dt = 0.04; ddt = 0.000000; T = 16.000; xn-xf = 0.0000; improve = 0.0000001; --> U^2*dt=2.6934;
% U = 7.3598 ; dt = 0.05; ddt = 0.000000; T = 20.000; xn-xf = 0.1470; improve = 0.0823429; --> U^2*dt=2.7083;
% U = 6.4222 ; dt = 0.066225; ddt = 0.000005; T = 26.490; xn-xf = 0.0050; improve = 0.0031986; --> U^2*dt=2.7314;
% 16.1308*0.013 = 3.3826;
% 14.6252*0.014 = 2.9946;
% 13.7285*0.015 = 2.8271;
% 13.1690*0.016 = 2.7748;
% 
% U = 14.1984; dt = 0.0137;T = 6.864;--> U^2*dt = 2.7673;
% Adjusting dt:
% U = 19.2459; dt = 0.01 ; T = 5;   --> U^2*dt = 3.7041;
% U = 16.6949; dt = 0.011; T = 5.5; --> U^2*dt = 3.0659;
% U = 15.3489; dt = 0.012; T = 6;   --> U^2*dt = 2.8271;
% U = 14.5991; dt = 0.013; T = 6.5; --> U^2*dt = 2.7707;
% U = 14.3176; dt = 0.0135;T = 6.75;--> U^2*dt = 2.7674;
% U = 14.2647; dt = 0.0136;T = 6.8; --> U^2*dt = 2.7673 394 717 2882;
% U = 14.2125; dt = 0.0137;T = 6.85;--> U^2*dt = 2.7673 203 713 3864;
% U = 14.2125; dt = 0.0137;T = 6.86;--> U^2*dt = 2.7673 199 886 3416;
% U = 14.2125; dt = 0.0137;T = 6.87;--> U^2*dt = 2.7673 199 889 3763; <---
% U = 14.1609; dt = 0.0138;T = 6.9; --> U^2*dt = 2.7673 192 608 1397;
% U = 14.0591; dt = 0.014; T = 7;   --> U^2*dt = 2.7672;
% U = 13.5650; dt = 0.015; T = 7.5; --> U^2*dt = 2.7601;
% U = 13.0897; dt = 0.016; T = 8;   --> U^2*dt = 2.7414;
% U = 12.6502; dt = 0.017; T = 8.5; --> U^2*dt = 2.7205;
% U = 12.2595; dt = 0.018; T = 9;   --> U^2*dt = 2.7053;
% U = 11.9159; dt = 0.019; T = 9.5; --> U^2*dt = 2.6978;
% U = 11.6093; dt = 0.02 ; T = 10;  --> U^2*dt = 2.6955;
% U = 9.4751;  dt = 0.03 ; T = 15;  --> U^2*dt = 2.6933;
% U = 8.2241;  dt = 0.04 ; T = 20;  --> U^2*dt = 2.7054;



x_target = [0;0];
u = zeros(iter_max,1); 
while true
    x = x0;
    x_traj = [];

    for iter = 1:iter_max
        % Prepare A,B to later calculate H
        R_big = JPhi(dt,x(1),x(2),u(iter));        
        A_store{iter} = R_big(1:2,1:2);       % A in Ax+Bu 
        B_store{iter} = R_big(1:2,3:3);       % B in Ax+Bu

        % The actual trajectory using adaptive step size mechanism
        [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 dt],[x;u(iter)]); 
        x = x_trajJ_fine(end,:)'; % the end of this sequence is x[k+1]
        x = x(1:2);
        x_traj = [x_traj x];      
    end
    
    error = norm(x-x_target);
    if error <= 0.0001
        break
    end
    % Calculate H
    H = B_store{1};
    for iter = 2:iter_max
        H = [A_store{iter}*H, B_store{iter}];
    end
    u = u - (H'*H+0.01*error^2*eye(iter_max))\(H'*(x-x_target));
end 
display('done1')
u1 = u;

% Part 2:
% Fig = openfig('vanderpol.fig');
% axis equal
% hold on 
% grid on
% step2_current_traj = scatter([x0(1),x_traj(1,:)],[x0(2),x_traj(2,:)],10,'k','filled');

u = u1;
%%
dt_max = 0.026;
mu = 50;
old_cost = norm(u)^2 + 1;

for iter2 = 1:10000      
    x = x0;    
    x_traj = [];
    for iter = 1:iter_max
        % Prepare A,B to later calculate H
        R_big = JPhi(dt,x(1),x(2),u(iter));  
        S_big = dPhi_dh_syms(dt,x(1),x(2),u(iter));
        % A, B, C in  Adx + Bdu + Cdt 
        A_store{iter} = R_big(1:2,1:2);       
        B_store{iter} = R_big(1:2,3:3);      
        C_store{iter} = S_big(1:2);
        % The actual trajectory using adaptive step size mechanism
        [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 dt],[x;u(iter)]); 
        x = x_trajJ_fine(end,:)'; % the end of this sequence is x[k+1]
        x = x(1:2);
        x_traj = [x_traj x];         
    end
    
    improve = old_cost - norm(u)^2;
    if iter2>1
        data = [norm(u), dt, du_and_ddt(end), dt*iter_max, norm(x-x_target), norm(u)^2*dt, improve, iter2];
        formatSpec = 'U = %.4f; dt = %.4f; ddt = %.6f; T = %.3f; xn-xf = %.4f; cost = %.4f; improve = %.7f; --> %.0f;\n';
        fprintf(formatSpec, data)
    end
    if improve < 1e-3
        mu=10;
    end
    if  norm(x-x_target) > 0.0050
        mu=100;
    end
    if  norm(x-x_target) > 0.0070
        mu=1000;
    end
    if improve < 1e-7 || norm(x-x_target)>0.01
        break
    end
    old_cost = norm(u)^2;
    
%     delete(step2_current_traj);
% %     step2_current_traj = plot([x0(1),x_traj(1,:)],[x0(2),x_traj(2,:)],'k','LineWidth',2);
%     step2_current_traj = scatter([x0(1),x_traj(1,:)],[x0(2),x_traj(2,:)],10,'r','filled');
%     drawnow 

    % create H:
    Hu = B_store{1};
    Ht = C_store{1};
    for iter = 2:iter_max
        Hu = [A_store{iter}*Hu, B_store{iter}];  
        Ht = [A_store{iter}*Ht, C_store{iter}]; 
    end
    H = [Hu,Ht*ones(iter_max,1)];

    options = optimset('display','off');
    du_and_ddt = quadprog( blkdiag((1+mu)*eye(iter_max),0), 2*[u',0], [zeros(2,iter_max),[1;-1]], [dt_max-dt;dt], H, x_target-x,[],[],[],options);
    u = u + du_and_ddt(1:iter_max);
    dt = dt + du_and_ddt(end);
end
display('done2') 
norm(u)^2*dt
% u_optimal = u;
% T_optimal = T;
% trajectory_optimal = x_traj;
% toc 

%%
Fig = openfig('vanderpol_inner2.fig');
axis equal
hold on 

dt = 0.005;
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
plot_discretized_trajectory = scatter(x_traj(1,:),x_traj(2,:),6, [0.2 0.1 0],'filled');
% plot_trajectory_optimal = scatter(trajectory_optimal(1,:),trajectory_optimal(2,:),15, 'r','filled');



function [Phi,Psi_p,JPhi,dPhi_dh_syms] = compute_Phi_and_JPhi2(p,F,xi,h)
Id = [];

for i = 1:length(xi)
    Id = [Id; xi(i)];
end

Psi = Id;             % [xi1; xi2; xi3; xi4; xi5] 
J = eye(length(xi));  % Identity 

Phi_syms = Psi;       % [xi1; xi2; xi3; xi4; xi5] 
JPhi_syms = J;        % Identity 
dPhi_dh_syms = 0*Psi;

for k_deg = 1:p
    tic
    Psi = J*F;
    Psi = simplify(Psi); % Calculate Psi from eq(1)
    
    J = jacobian(Psi,xi);
    J = simplify(J);   % Jacobian of Psi
    
    % Taylor approxi. of the flow (discretization purpose)
    Phi_syms = Phi_syms + (h^k_deg/factorial(k_deg))*Psi;
    % Taylor approxi. of the jacobian of flow (linearization purpose)
    JPhi_syms = JPhi_syms + (h^k_deg/factorial(k_deg))*J;
    % Partial flow/ partial h
    dPhi_dh_syms = dPhi_dh_syms + (h^(k_deg-1)/factorial(k_deg-1))*Psi;
    toc
end
% Turn symbolic expressions into functions 
Phi = matlabFunction(simplify(Phi_syms),'Vars',[h xi]);   % The flow 
Psi_p = matlabFunction(simplify(Psi),'Vars',xi);          % The last Psi 
JPhi = matlabFunction(simplify(JPhi_syms),'Vars',[h xi]); % Jacobian of the flow
dPhi_dh_syms = matlabFunction(simplify(dPhi_dh_syms),'Vars',[h xi]); % Temporal sensitivity of the flow
end



