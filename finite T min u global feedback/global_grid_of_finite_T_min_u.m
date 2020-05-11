clear
clc
close all

syms xi1 xi2 xi3 dt

F = [xi2;
     -xi1+xi2-xi1^2*xi2+xi3;
      0];
  
F = simplify(F);
     
p = 6;
tic
[Phi,Psi_p,JPhi] = compute_Phi_and_JPhi(p,F,[xi1 xi2 xi3],dt);
toc
%%
% x0 = [-2.4;3]; T=5
x0 = [-3;3]; %T=3  
% x0 = [-1.25189390985363;1.54007819071531];
x_target = [0;0];
T = 4;    %5-64*0.01; %3-36*0.01;
dt = 0.01;
iter_max = ceil(T/dt);
u = zeros(iter_max,1); 
while true
    x = x0;
    x_traj = [];
    for iter = 1:iter_max
        % Prepare A,B to later calculate H
        R_big = JPhi(dt,x(1),x(2),u(iter));        
        R_store{iter} = R_big(1:2,1:2);       % A in Ax+Bu 
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
        H = [R_store{iter}*H, B_store{iter}];
    end
    u = u - (H'*H+0.001*error^2*eye(iter_max))\(H'*(x-x_target));
end 

old_cost = norm(u);
mu = 5;
for iter2 = 1:500       
    x = x0;    
    x_traj = [];
    for iter = 1:iter_max
        % Prepare A,B to later calculate H
        R_big = JPhi(dt,x(1),x(2),u(iter));        
        R_store{iter} = R_big(1:2,1:2);       % A in Ax+Bu 
        B_store{iter} = R_big(1:2,3:3);       % B in Ax+Bu

        % The actual trajectory using adaptive step size mechanism
        [~, x_trajJ_fine] = adaptive_taylor(p,Phi,Psi_p,[0 dt],[x;u(iter)]); 
        x = x_trajJ_fine(end,:)'; % the end of this sequence is x[k+1]
        x = x(1:2);
        x_traj = [x_traj x];         
    end
    [norm(x_target-x) norm(u)]
    % create H:
    H = B_store{1};
    for iter = 2:iter_max
        H = [R_store{iter}*H, B_store{iter}];       
    end
    options = optimset('display','off');
    u = u + quadprog((1+mu)*eye(iter_max),2*u,[],[],H,x_target-x,[],[],[],options);
    improve = old_cost - norm(u);
    old_cost = norm(u);
    if improve<1e-4
        break 
    end
end
data = [x0(1) x0(2) T norm(u) norm(x_target-x) iter2];
formatSpec = '(%.2f;%.2f); T = %.2f; U = %.4f; xn-xf = %.4f;  --> %.0f;\n';
fprintf(formatSpec, data)    

Fig = openfig('vanderpol.fig');
axis equal
hold on 
grid on
step2_current_traj = scatter([x0(1),x_traj(1,:)],[x0(2),x_traj(2,:)],10,'r','filled');
% delete(step2_current_traj)
