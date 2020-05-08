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

%% Part 1:
tic
dt = 0.01;
T = 4;
iter_max = ceil(T/dt);


x0 = [-2;3];%U = 26.6083; t = 0.01; dt/t = -0.000; T = 4.075; cost = 23.82; xn-xf = 0.0000; improve = 0.0000001; --> 624;

x0 = [1.978;1.901]; 

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
axis equal
hold on 
grid on
step2_current_traj = scatter([x0(1),x_traj(1,:)],[x0(2),x_traj(2,:)],10,'k','filled');

u = u1;
t = dt; % t is the time interval for each step 
t_max = 0.03;
mu = 1;
old_cost = t*norm(u)^2 + (iter_max*t)^2 + 1;

for iter2 = 1:10000
       
    x = x0;    
    x_traj = [];
    for iter = 1:iter_max
        % Prepare A,B to later calculate H
        R_big = JPhi(t,x(1),x(2),u(iter));  
        S_big = dPhi_dh_syms(t,x(1),x(2),u(iter));
        % A, B, C in  Adx + Bdu + Cdt 
        A_store{iter} = R_big(1:2,1:2);       
        B_store{iter} = R_big(1:2,3:3);      
        C_store{iter} = S_big(1:2);

        % The actual trajectory using adaptive step size mechanism
        [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 t],[x;u(iter)]); 
        x = x_trajJ_fine(end,:)'; % the end of this sequence is x[k+1]
        x = x(1:2);
        x_traj = [x_traj x];         
    end
    
    improve = old_cost - (t*norm(u)^2 + (iter_max*t)^2);
    if iter2>1
        data = [norm(u),t, 100*du_and_dt(end)/t, t*iter_max, old_cost, norm(x-x_target), improve, iter2];
        formatSpec = 'U = %.4f; t = %.2f; dt/t = %.3f; T = %.3f; cost = %.2f; xn-xf = %.4f; improve = %.7f; --> %.0f;\n';
        fprintf(formatSpec, data)
    end
    
    if  norm(x-x_target) > 0.0050
        mu=100;
    end
    if improve < 1e-3
        mu=1;
    end
    if improve < 1e-7 || norm(x-x_target)>0.01
        break
    end
    old_cost = t*norm(u)^2 + (iter_max*t)^2;
    
    delete(step2_current_traj);
%     step2_current_traj = plot([x0(1),x_traj(1,:)],[x0(2),x_traj(2,:)],'k','LineWidth',2);
    step2_current_traj = scatter([x0(1),x_traj(1,:)],[x0(2),x_traj(2,:)],10,'k','filled');
    drawnow 

    % create H:
    Hu = B_store{1};
    Ht = C_store{1};
    for iter = 2:iter_max
        Hu = [A_store{iter}*Hu, B_store{iter}];  
        Ht = [A_store{iter}*Ht, C_store{iter}]; 
    end
    H = [Hu,Ht*ones(iter_max,1)];

    options = optimset('display','off');
    du_and_dt = quadprog( blkdiag((t+mu)*eye(iter_max),iter_max^2), 2*[t*u',iter_max^2*t], [zeros(2,iter_max),[1;-1]], [t_max-t;t], H, x_target-x,[],[],[],options);
    u = u + du_and_dt(1:iter_max);
    t = t + du_and_dt(end);
end
display('done2')
u_optimal = u;
T_optimal = T;
trajectory_optimal = x_traj;
toc 

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



