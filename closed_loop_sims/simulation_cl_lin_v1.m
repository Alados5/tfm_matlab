%{
Closed-loop simulation of an MPC applied to a cloth model
Original code by David Parent
Modified by Adria Luque so both COM and SOM are linear models,
as a first step to adapt the code into C++ to apply it to the real robot
%}
clear; close all; clc;

addpath('..\required_files\cloth_model_DParent')
addpath('..\required_files\casadi-toolbox')
import casadi.*

plotAnim = 0;

% General Parameters
nCOM = 4;
nSOM = 4;
NTraj = 2;
Ns_og = 10;
Hp = 30;
Ts = 0.010;

% Opti parameters
xbound = 1;
ubound = 1*1e-3;
gbound = 0;  % 0 -> Equality constraint
W_Q = 1;
W_T = 1;
W_R = 10;

% -------------------

% Define COM parameters
nxC = nCOM;
nyC = nCOM;
COMlength = nxC*nyC;
COM = struct;
COM.row = nxC;
COM.col = nyC;
COM.mass = 0.1;
COM.grav = 9.8;
COM.dt = Ts;

% Parameters from calibration
%{
% For a 4x4 (Optimized)
COM.stiffness = [-305.6028  -13.4221 -225.8987]; 
COM.damping = [-4.0042   -2.5735   -3.9090]; 
COM.z_sum = 0.0312;
%}
%{
% For a 4x4 (Learnt Exp2)
COM.stiffness = [-302.1494 -13.7817 -220.1087]; 
COM.damping = [-4.2148 -1.0375 -3.6974]; 
COM.z_sum = 0.0291;
%}
%{
% For a 4x4 (Learnt Exp4)
COM.stiffness=[-247.4008 -201.2089 -405.4698];
COM.damping=[-3.4821   -3.1614   -4.5056];
COM.z_sum=0.0150;
%}
%{x
% For a 4x4 (Learnt Exp7)
COM.stiffness=[-331.6361 -27.5886 -442.8003];
COM.damping=[-4.2262 -2.6943 -4.8822];
COM.z_sum=0.0142;
%}

% Controlled coordinates (upper corners in x,y,z)
COM_node_ctrl = [nxC*(nyC-1)+1, nxC*nyC];
COM.coord_ctrl = [COM_node_ctrl, ...
                  COM_node_ctrl+nxC*nyC, ...
                  COM_node_ctrl+2*nxC*nyC];

% Define the SOM (LINEAR)
nxS = nSOM;
nyS = nSOM;
SOMlength = nxS*nyS;
SOM = struct;
SOM.row = nxS;
SOM.col = nyS;
SOM.mass = 0.1;
SOM.grav = 9.8;
SOM.dt = Ts;
pos = load(['../data/pos', num2str(Ns_og), '_', num2str(NTraj), '.mat']);
pos = pos.pos;

% Reduce initial position if SOM size is not 10
i = 0:nSOM^2-1;
redc = mod(i,nSOM)*(Ns_og-1)/(nSOM-1) + ...
       floor(i/nSOM)*(Ns_og-1)/(nSOM-1)*Ns_og + 1;
posr = pos(redc,:);


% Parameters from calibration
%{x
% For a 4x4 (Optimized)
SOM.stiffness = [-305.6028  -13.4221 -225.8987];%.*[1 0.8 2.0]; 
SOM.damping = [-4.0042   -2.5735   -3.9090];%.*[1 0.8 1.25]; 
SOM.z_sum = 0.0312;%*1.5;
%}
%{
% For a 4x4 (Learnt)
SOM.stiffness = [-302.1494 -13.7817 -220.1087]; 
SOM.damping = [-4.2148 -1.0375 -3.6974]; 
SOM.z_sum = 0.0291;
%}
%{
% For a 4x4 (Learnt Exp4)
SOM.stiffness=[-247.4008 -201.2089 -405.4698];
SOM.damping=[-3.4821   -3.1614   -4.5056];
SOM.z_sum=0.0150;
%}
%{
% For a 10x10
SOM.stiffness = [-300  -10.0001 -200.0021]; 
SOM.damping = [-4.0001   -2.0002   -3.7444]; 
SOM.z_sum = 0.3471;
%}
%{
% For a 7x7
SOM.stiffness = [-535.9092    2.7121 -289.3164]; 
SOM.damping = [-5.2901   -1.2518   -4.0844]; 
SOM.z_sum = 0.0713;
%}

% Controlled coordinates (upper corners in x,y,z)
SOM_node_ctrl = [nxS*(nyS-1)+1, nxS*nyS];
SOM.coord_ctrl = [SOM_node_ctrl SOM_node_ctrl+nxS*nyS SOM_node_ctrl+2*nxS*nyS];

% Define initial position of the nodes (needed for ext_force)
% Second half is velocity (initial v=0)
x_ini_SOM = [reshape(posr,[3*nxS*nyS 1]); zeros(3*nxS*nyS,1)];
[reduced_pos,~] = take_reduced_mesh(x_ini_SOM(1:3*nxS*nyS),x_ini_SOM(3*nxS*nyS+1:6*nxS*nyS), nSOM, nCOM);
x_ini_COM = [reduced_pos; zeros(3*nxC*nyC,1)];

% Initial position of the nodes
SOM.nodeInitial = lift_z(reshape(x_ini_SOM(1:3*nxS*nyS), [nxS*nyS,3]), SOM);
COM.nodeInitial = lift_z(reshape(x_ini_COM(1:3*nxC*nyC), [nxC*nyC,3]), COM);


% Find initial spring length in each direction x,y,z
[SOM.mat_x, SOM.mat_y, SOM.mat_z] = compute_l0_linear(SOM,0);
[COM.mat_x, COM.mat_y, COM.mat_z] = compute_l0_linear(COM,0);

% Find linear matrices
[A_SOM, B_SOM, f_SOM] = create_model_linear_matrices(SOM);
[A_COM, B_COM, f_COM] = create_model_linear_matrices(COM);

%% Start casADi optimization problem
% Declare model variables
phi = SX.sym('phi',3*COM.row*COM.col);
dphi = SX.sym('dphi',3*COM.row*COM.col);
x = [phi; dphi];
n_states = length(x); %should be 96 in a 4x4 model
u = SX.sym('u',6);

% Define model equations
x_next = SX.sym('xdot',6*COM.row*COM.col);
x_next(:) = A_COM*x + B_COM*u + COM.dt*f_COM;

% (x,u)->(x_next)
stfun = Function('stfun',{x,u},{x_next}); % nonlinear mapping function f(x,u)

% Lower corner coordinates for both models
coord_l = [1 nyC 1+nxC*nyC nxC*nyC+nyC 2*nxC*nyC+1 2*nxC*nyC+nyC]; 
coord_nl = [1 nxS 1+nxS*nyS nxS*nyS+nxS 2*nxS*nyS+1 2*nxS*nyS+nxS];

% Parameters in the optimization problem: initial state and reference
P = SX.sym('P',n_states,Hp+1); 

% Start with an empty NLP
w=[]; %variables to optimize
lbw = []; %lower bounds
ubw = []; %upper bounds
obj = 0; %objective function
g=[]; %constraints
lbg = [];
ubg = [];

% Initial condition (initial state)
Xk = P(:,1);

% Weigths calculation (adaptive): direction pointing from the actual
% lower corner position to the desired position at the end of the interval
ln = P([1 3 5],end) - Xk(coord_l([1,3,5])); %ln = left node
ln = abs(ln)./(norm(ln)+10^-6);
rn = P([2,4,6],end) - Xk(coord_l([2,4,6])); %rn = right node
rn = abs(rn)./(norm(rn)+10^-6);
Q = diag([ln(1) rn(1) ln(2) rn(2) ln(3) rn(3)]);

take_x=[];
take_u=[];
i=1;

for k = 1:Hp
    % Artificial reference
    r_a = SX.sym(['r_a_' num2str(k)],6);
    w = [w; r_a];
    
    lbw = [lbw; -xbound*ones(6,1)]; 
    ubw = [ubw;  xbound*ones(6,1)];
    
    take_x=[take_x;(i:(i+6-1))'];
    i=i+6;
    
    % New NLP variable for the control
    Uk = SX.sym(['U_' num2str(k)],6);
    w = [w; Uk];
    % The bounds of the control action are very relevant! 
    %   If too large, the COM allows it and the SOM will not
    %   (the COM is a spring which can be infinetely stretched in a step)
    lbw = [lbw; -ubound*ones(6,1)]; 
    ubw = [ubw;  ubound*ones(6,1)];

    take_u=[take_u;(i:(i+6-1))'];
    i=i+6;
    
    % Integrate till the end of the interval
    Xk_next = stfun(Xk,Uk);
    Xkn_r = Xk_next(coord_l);
    
    % Update objective function
    obj = obj + W_Q*(Xkn_r-r_a)'*Q*(Xkn_r-r_a);
    obj = obj + W_T*(P(1:6,k+1)-r_a)'*(P(1:6,k+1)-r_a);
    obj = obj + W_R*(Uk'*Uk);

    Xk = Xk_next;
end

% Terminal constraint
g = [g; Xkn_r-r_a ]; 
lbg = [lbg; -gbound*ones(6,1)];
ubg = [ubg;  gbound*ones(6,1)];

nlp_prob = struct('f', obj, 'x', w, 'g', g, 'p', P);
opts = struct;
opts.print_time = 0;
opts.ipopt.print_level = 0; %0 to print the minimum, 3 to print the maximum
opts.ipopt.warm_start_init_point = 'yes'; %warm start

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP

%%
% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------

% Define trajectory to follow
phi_l_Traj = load(['../data/trajectories/phi_',num2str(NTraj),'L.csv']);
phi_r_Traj = load(['../data/trajectories/phi_',num2str(NTraj),'R.csv']);

u_ini = x_ini_SOM(SOM.coord_ctrl);

% Initialize things
reference = zeros(6*COM.row*COM.col, Hp+1);
store_state(:,1) = x_ini_SOM;
store_u(:,1) = zeros(6,1);

t0 = tic;
printX = 10;
for tk=2:size(phi_l_Traj,1)
    
    % The last Hp timesteps, trajectory should remain constant
    if tk>=size(phi_l_Traj,1)-Hp 
        Traj_l_Hp = repmat(phi_l_Traj(end,:), Hp,1);
        Traj_r_Hp = repmat(phi_r_Traj(end,:), Hp,1);
    else
        Traj_l_Hp = phi_l_Traj(tk:tk+Hp-1,:);
        Traj_r_Hp = phi_r_Traj(tk:tk+Hp-1,:);
    end
    
    % Define reference in the prediction horizon (sliding window)
    reference(:,1) = x_ini_COM;
    reference(1,2:end) = Traj_l_Hp(:,1)';
    reference(2,2:end) = Traj_r_Hp(:,1)';
    reference(3,2:end) = Traj_l_Hp(:,2)';
    reference(4,2:end) = Traj_r_Hp(:,2)';
    reference(5,2:end) = Traj_l_Hp(:,3)';
    reference(6,2:end) = Traj_r_Hp(:,3)';
    
    % Initial seed of the optimization (for "u" and "r^a")
    args_x0  = repmat([reference(1:6,end)+[0;0;0;0;0.3;0.3];zeros(6,1)],Hp,1);
    
    % Find the solution "sol"
    sol = solver('x0', args_x0, 'lbx', lbw, 'ubx', ubw, ...
                 'lbg', lbg, 'ubg', ubg, 'p', reference);

    % Get only controls from the solution
    u = reshape(full(sol.x(take_u)),6,Hp)'; 
    u_lin = u(1,1:end)';
    
    next_state_SOM = A_SOM*store_state(:,tk-1) + B_SOM*u_lin + SOM.dt*f_SOM;
    
    phi_ini_SOM = next_state_SOM(1:3*nxS*nyS);
    dphi_ini_SOM = next_state_SOM((1+3*nxS*nyS):6*nxS*nyS);    
        
    % Close the loop
    [phired, dphired] = take_reduced_mesh(phi_ini_SOM,dphi_ini_SOM, nSOM, nCOM);
    x_ini_COM = [phired; dphired];
    
    % Store things
    store_state(:,tk) = [phi_ini_SOM; dphi_ini_SOM];
    store_u(:,tk) = u_lin;
    
    if(mod(tk,printX)==0)
        t10 = toc(t0)*1000;
        fprintf(['Iter: ', num2str(tk), ...
            ' \t Avg. time/iter: ', num2str(t10/printX), ' ms \n']);
        t0 = tic;
    end
end

%% PLOT LOWER CORNERS
time = 0:Ts:size(store_state,2)*Ts-Ts;

fig1 = figure(1);
fig1.Color = [1,1,1];

subplot(1,2,1);
plot(time,store_state(coord_nl([1 3 5]),:)', 'linewidth',1.5)
hold on
plot(time, phi_l_Traj, '--k', 'linewidth',1.2)
hold off
title('\textbf{Left lower corner}', 'Interpreter', 'latex')
legend('x','y','z','Ref', 'Interpreter', 'latex')
grid on
xlabel('t [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(1,2,2)
plot(time,store_state(coord_nl([2 4 6]),:)', 'linewidth',1.5)
hold on
plot(time, phi_r_Traj, '--k', 'linewidth',1.2)
hold off
title('\textbf{Right lower corner}', 'Interpreter', 'latex')
legend('x','y','z','Ref', 'Interpreter', 'latex')
grid on
xlabel('t [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

avg_error_1 = 1000*[mean(abs(store_state(coord_nl(1),:)'-phi_l_Traj(1:end,1)))
                    mean(abs(store_state(coord_nl(3),:)'-phi_l_Traj(1:end,2)))
                    mean(abs(store_state(coord_nl(5),:)'-phi_l_Traj(1:end,3)))];
avg_error_4 = 1000*[mean(abs(store_state(coord_nl(2),:)'-phi_r_Traj(1:end,1)))
                    mean(abs(store_state(coord_nl(4),:)'-phi_r_Traj(1:end,2)))
                    mean(abs(store_state(coord_nl(6),:)'-phi_r_Traj(1:end,3)))];
avg_error = mean([avg_error_1 avg_error_4]);

display(avg_error_1);
display(avg_error_4);
display(avg_error);



%% PLOT UPPER CORNERS
fig2 = figure(2);
fig2.Color = [1,1,1];
fig2.Position = fig2.Position + [100 0 0 0];

subplot(1,2,1);
plot(time,store_state(SOM.coord_ctrl([1 3 5]),:)','linewidth',1.5)
title('\textbf{Left upper corner}', 'Interpreter', 'latex')
legend('x','y','z', 'Interpreter', 'latex')
grid on
xlabel('t [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(1,2,2)
plot(time,store_state(SOM.coord_ctrl([2 4 6]),:)','linewidth',1.5)
title('\textbf{Right upper corner}', 'Interpreter', 'latex')
legend('x','y','z', 'Interpreter', 'latex')
grid on
xlabel('t [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');


%% PLOT CLOTH MOVING
fig3 = figure(3);
fig3.Color = [1,1,1];
fig3.Position = fig3.Position + [200 -60 60 60];

store_pos = store_state(1:3*SOMlength,:);

store_x = store_pos(1:SOMlength,:);
limx = [floor(min(store_x(:))*10), ceil(max(store_x(:))*10)]/10;
store_y = store_pos(SOMlength+1:2*SOMlength,:);
limy = [floor(min(store_y(:))*10), ceil(max(store_y(:))*10)]/10;
store_z = store_pos(2*SOMlength+1:3*SOMlength,:);
limz = [floor(min(store_z(:))*10), ceil(max(store_z(:))*10)]/10;

scatter3(store_x(:,1), store_y(:,1), store_z(:,1), '.b');
hold on
scatter3(store_x(1,1), store_y(1,1), store_z(1,1), '*b');
hold off
axis equal;
xlim(limx);
ylim(limy);
zlim(limz);
set(gca, 'TickLabelInterpreter','latex');
xlabel('x', 'Interpreter','latex');
ylabel('y', 'Interpreter','latex');
zlabel('z', 'Interpreter','latex');

if(plotAnim)
    %hold on;
    for t=2:size(store_state,2)

        scatter3(store_x(:,t), store_y(:,t), store_z(:,t), '.b');
        axis equal
        xlim(limx);
        ylim(limy);
        zlim(limz);
        set(gca, 'TickLabelInterpreter','latex');
        xlabel('x', 'Interpreter','latex');
        ylabel('y', 'Interpreter','latex');
        zlabel('z', 'Interpreter','latex');

        pause(1e-6);
    end
    hold off
end
