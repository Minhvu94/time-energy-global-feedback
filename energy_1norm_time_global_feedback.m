clear
clc
close all
clear all 
syms xi1 xi2 xi3 dt

F = [xi2;
     -xi1+xi2-xi1^2*xi2+xi3;
      0];
  
F = simplify(F);
     
p = 6;
[Phi,Psi_p,JPhi,dPhi_dh_syms] = compute_Phi_and_JPhi2(p,F,[xi1 xi2 xi3],dt);

%% Part 1:
dt = 0.01;
T = 2;
iter_max = ceil(T/dt);

% x0 = [2;-3]; % for vanderpol1: [-0.00772105272954248;-3.15798353408161]; [-2.22589498649188;-0.623240650204037]; [-1.68115398300972;0.720565630180224]
% x0 = [2;3]; % for vanderpol2: [2.41258997419395;-0.391530160531988];  [2.03887555867469;-0.579682221532998];  [1.79955404232931;-0.673143387123778]
x0 = [2;3];

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
Fig = openfig('vanderpol2.fig');
axis equal
hold on 
grid on
plot(x0(1),x0(2),'xr','LineWidth',2);
step2_current_traj = plot([x0(1),x_traj(1,:)],[x0(2),x_traj(2,:)],'k','LineWidth',2);
u = u1;
T = dt*ones(iter_max,1);
%%
mu = 100;

for iter2 = 1:3000
       
    x = x0;    
    x_traj = [];
    for iter = 1:iter_max
        % Prepare A,B to later calculate H
        R_big = JPhi(T(iter),x(1),x(2),u(iter));  
        S_big = dPhi_dh_syms(T(iter),x(1),x(2),u(iter));
        % A, B, C in  Adx + Bdu + Cdt 
        A_store{iter} = R_big(1:2,1:2);       
        B_store{iter} = R_big(1:2,3:3);      
        C_store{iter} = S_big(1:2);

        % The actual trajectory using adaptive step size mechanism
        [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 T(iter)],[x;u(iter)]); 
        x = x_trajJ_fine(end,:)'; % the end of this sequence is x[k+1]
        x = x(1:2);
        x_traj = [x_traj x];         
    end

    data = [norm(u), sum(T), norm(x-x_target), iter2];
    formatSpec = 'U = %.4f; T = %.3f; cost = %.4f; --> %.0f;\n';
    fprintf(formatSpec, data)
    
%     figure(1)
    delete(step2_current_traj);
    step2_current_traj = plot([x0(1),x_traj(1,:)],[x0(2),x_traj(2,:)],'k','LineWidth',2);
%     figure(2)
%     hold off
%     plot(T,'b','LineWidth',2);
%     grid on
    drawnow 

    % create H:
    Hu = B_store{1};
    Ht = C_store{1};
    for iter = 2:iter_max
        Hu = [A_store{iter}*Hu, B_store{iter}];  
        Ht = [A_store{iter}*Ht, C_store{iter}]; 
    end
    H = [Hu,Ht];

    options = optimset('display','off');
    du_and_dt = quadprog( blkdiag((0.001+mu)*eye(iter_max),mu*eye(iter_max)), [0.001*u',ones(1,iter_max)], [zeros(iter_max),-1*eye(iter_max)], T, H, x_target-x,[],[],[],options);
    u = u + du_and_dt(1:iter_max);
    T = T + du_and_dt(iter_max+1:end);
    
%     pause(1)

end

display('done2')
% data = [data [u(1);T - decrease]];
u_optimal = u;

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



