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

%%
[X,Y] = meshgrid(-3:1:3 , -3:1:3);
grids = [X(:)  Y(:)];

for rep=1:length(grids)
tic
fprintf('rep %.0f: \n',rep)
x0 = grids(rep,:)';
x_target = [0;0];
if x0(1) == 0 && x0(2) == 0
    display('Skip the origin!')
else   
    %  Part 1:
    dt = 0.01;
    T = 3;
    iter_max = ceil(T/dt);
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
        u = u - (H'*H+0.001*error^2*eye(iter_max))\(H'*(x-x_target));
    end 
    display('done 1')
    u1 = u;

    % Part 2:
    u = u1;
    T = dt*ones(iter_max,1);
    alpha = 1;
    beta = 1;
    mu = 100;
    old_cost = norm(u)^2+norm(T)^2+10;

    for iter2 = 1:7000
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

        improve = old_cost - (alpha*norm(u)^2+beta*norm(T)^2);
%         data = [norm(u), sum(T), old_cost, norm(x-x_target), improve, iter2];
%         formatSpec = 'U = %.4f; T = %.3f; U^2+T^2 = %.2f; xn-xf = %.4f; improve = %.6f; --> %.0f;\n';
%         fprintf(formatSpec, data)
        if improve<0.001
            mu=50;
        end
        if norm(x-x_target)>0.0050
            mu = 1000;
    %     elseif norm(x-x_target)<0.0020
    %         mu = 100;
        end
        if improve < 1e-5 || norm(x-x_target)>0.01
            break
        end
        old_cost = alpha*norm(u)^2+beta*norm(T)^2;

%         delete(step2_current_traj);
%         step2_current_traj = plot([x0(1),x_traj(1,:)],[x0(2),x_traj(2,:)],'k','LineWidth',2);
%         drawnow 

        % create H:
        Hu = B_store{1};
        Ht = C_store{1};
        for iter = 2:iter_max
            Hu = [A_store{iter}*Hu, B_store{iter}];  
            Ht = [A_store{iter}*Ht, C_store{iter}]; 
        end
        H = [Hu,Ht];

        options = optimset('display','off');
        du_and_dt = quadprog( blkdiag((alpha+mu)*eye(iter_max),(beta+mu)*eye(iter_max)), 2*[alpha*u',beta*T'], [zeros(iter_max),-1*eye(iter_max)], T, H, x_target-x,[],[],[],options);
        u = u + du_and_dt(1:iter_max);
        T = T + du_and_dt(iter_max+1:end);

    end

    fprintf('iter = %.0f; xn-xf = %.3f; U^2+T^2 = %.4f; improvement = %.6f; --> done2 \n',iter2,norm(x-x_target),old_cost,improve)

    store_u_optimal(:,rep) = u;
    store_T_optimal(:,rep) = T;
    store_trajectory_optimal{rep} = [x0,x_traj];
end
toc
end

save('u_optimal.mat','store_u_optimal')
save('T_optimal.mat','store_T_optimal')
save('trajectory_optimal.mat','store_trajectory_optimal')

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



