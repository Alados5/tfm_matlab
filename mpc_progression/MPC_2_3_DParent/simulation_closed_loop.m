% Solve the MPC problem where the COM is a linear cloth model and the SOM
% is the Coltraro's cloth model.
%
% Author: David Parent, davidparentalonso@gmail.com
% Last review: 07/02/2021

clear; close all; clc;

addpath('required_files\cloth_model_FColtraro')
addpath('required_files\cloth_model_DParent')
addpath('..\..\required_files\casadi-toolbox')
import casadi.*

% Define COM parameters(computation oriented model)
COM = struct;
COM.row = 4; COM.col = 4;
COM.mass = 0.1;
COM.grav = 9.8;
COM.dt = 0.01;

% parameters from calibration
COM.stiffness = [-305.6028  -13.4221 -225.8987]; 
COM.damping = [-4.0042   -2.5735   -3.9090]; 
COM.z_sum = 0.0312;

% Define the SOM
nx = 10; ny = nx;
[SOM,pos] = initialize_model_parameters(nx,ny);

% Definite initial position of the nodes (defined here because we need it
% to compute ext_force
x_ini_SOM = [reshape(pos,[3*nx*ny 1]);zeros(3*nx*ny,1)]; %initial velocity=0
[reduced_pos,~] = take_reduced_mesh(x_ini_SOM(1:300),x_ini_SOM(301:600));
x_ini_COM = [reduced_pos;zeros(48,1)];

% Initial position of the COM nodes
COM.nodeInitial = lift_z(reshape(x_ini_COM(1:48), [16,3]), COM);

% Find initial spring length in each direction x,y,z
[COM.mat_x, COM.mat_y, COM.mat_z] = compute_l0_linear(COM,0);

% Find linear matrices
[A, B, ext_force] = create_model_linear_matrices(COM);

%% Start casADi optimization problem

% Declare model variables
phi = SX.sym('phi',3*COM.row*COM.col);
dphi = SX.sym('dphi',3*COM.row*COM.col);
x = [phi; dphi];
u = SX.sym('u',6);
n_states = length(x); %should be 96 in a 4x4 model

% Define model equations
%x_next = SX.sym('xdot',6*COM.row*COM.col);
%x_next(:) = ;

% (x,u)->(x_next)
f = Function('f',{x,u},{A*x + B*u + (COM.dt).*ext_force}); % nonlinear mapping function f(x,u)

nr = COM.row;
nc = COM.col;
cord_l = [1 nc 1+nr*nc nr*nc+nc 2*nr*nc+1 2*nr*nc+nc]; %coordenates of the linear model (of the lower corners)
cord_nl = [1 nx 1+nx*ny nx*ny+nx 2*nx*ny+1 2*nx*ny+nx]; %coordenates of the non-linear model (of the lower corners)
COM.coord_controlados = [13 16 29 32 45 48];

ubound = 0.8e-3; %0.8e-3, 50e-3
R = 10; %10, 20
N = 30; % Number of control intervals
P = SX.sym('P',n_states,N+1); %Parameters in the optimization problem. It includes initial state and reference

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
nl = P([1 3 5],end) - Xk(cord_l([1,3,5])); %nl == node left
nl = abs(nl)./(norm(nl)+10^-6);
nr = P([2,4,6],end) - Xk(cord_l([2,4,6])); %nr == node right
nr = abs(nr)./(norm(nr)+10^-6);
Q = diag([nl(1) nr(1) nl(2) nr(2) nl(3) nr(3)]);

take_x=[];
take_u=[];
i=1;

for k = 1:N
    % Artificial reference
    r_a = SX.sym(['r_a_' num2str(k)],6);
    w = [w; r_a];
    lbw = [lbw;(-1)*ones(6,1)];
    ubw = [ubw;(1)*ones(6,1)];
    take_x=[take_x;(i:(i+6-1))'];
    i=i+6;
    
    % New NLP variable for the control
    Uk = SX.sym(['U_' num2str(k)],6);
    w = [w; Uk];
    %The bounds of the control input are very relevant!! If too big, the COM
    %allows it but the SOM not (the COM is a spring that from a time step
    %to the next one can be mathematically infinetely stretched). The
    %bounds should be approximately the maximum displacement expected.
    lbw = [lbw; -ubound*ones(6,1)]; 
    ubw = [ubw;  ubound*ones(6,1)];

    take_u=[take_u;(i:(i+6-1))'];
    i=i+6;
    
    % Integrate till the end of the interval
    Xk_next = f(Xk,Uk);
    
    % Update objective function
    obj = obj + (1000*Xk_next(cord_l)- 1000*r_a)'*Q*(1000*Xk_next(cord_l)- 1000*r_a);
    obj = obj + R*1000*1000*(Uk'*Uk);
    obj = obj + (1000*P(1:6,k+1)-1000*r_a)'*(1000*P(1:6,k+1)-1000*r_a);

    Xk = Xk_next;
end

% Terminal constraint
g = [g; Xk_next(cord_l)-r_a ]; 
lbg = [lbg; zeros(6,1)]; %equality constraint
ubg = [ubg; zeros(6,1)];

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
phi_l_Traj = load('required_files\trajectories\phi_l_3D.csv');
phi_r_Traj = load('required_files\trajectories\phi_r_3D.csv');

u_ini = x_ini_SOM(SOM.coord_controlados);
u_bef = u_ini;

% Counters for the loop
sHp = 1; % start prediction horizon
Hp = N; % end prediction horizon

% Initialize things
reference = zeros(6*COM.row*COM.col,N+1);
store_state(:,1) = x_ini_SOM;
store_state(:,2) = x_ini_SOM;
store_u(:,1) = zeros(6,1);
store_u(:,2) = zeros(6,1);

tT = 0;
%t0 = tic;
for i=3:size(phi_l_Traj,1)
    Hp = Hp+1;
    sHp = sHp+1;

    if i>=size(phi_l_Traj,1)-N %the last N timesteps, trajectory should remain constant
        Traj_l_Hp = repmat(phi_l_Traj(end,:),N);
        Traj_r_Hp = repmat(phi_r_Traj(end,:),N);
    else
        Traj_l_Hp = phi_l_Traj(sHp:Hp,:);
        Traj_r_Hp = phi_r_Traj(sHp:Hp,:);
    end
    
    % Define reference in the prediction horizon (moving window)
    reference(:,1) = x_ini_COM;
    reference(1,2:end) = Traj_l_Hp(:,1)';
    reference(2,2:end) = Traj_r_Hp(:,1)';
    reference(3,2:end) = Traj_l_Hp(:,2)';
    reference(4,2:end) = Traj_r_Hp(:,2)';
    reference(5,2:end) = Traj_l_Hp(:,3)';
    reference(6,2:end) = Traj_r_Hp(:,3)';
    
    % Initial seed of the optimization (for "u" and "r^a")
    args.x0  = repmat([reference(1:6,end)+[0;0;0;0;0.3;0.3];zeros(6,1)],N,1);
    args.p = reference;
    
    % Find the solution "sol"
    t0 = tic;
    sol = solver('x0', args.x0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg, 'p', args.p);

    u = reshape(full(sol.x(take_u)),6,N)'; % get only controls from the solution 
    u_SOM = u(1,1:end)' + u_bef;
    u_bef = u_SOM; %update
    
    tT=tT+toc(t0);
    [next_state_SOM] = cloth_simulator_secondorder([store_state(:,i-1);store_state(:,i-2)],u_SOM,SOM);
    
    phi_ini_SOM = full(next_state_SOM(1:3*SOM.n_nodos));
    dphi_ini_SOM = full(next_state_SOM((1+3*SOM.n_nodos):6*SOM.n_nodos));    
    %t0 = tic;
    
    % Close the loop
    [phired,dphired] = take_reduced_mesh(phi_ini_SOM,dphi_ini_SOM);
    x_ini_COM = [phired;dphired];
    
    % Store things
    store_state(:,i) = [phi_ini_SOM;dphi_ini_SOM];
    store_u(:,i) = u(1,1:end)';    
    
    if(mod(i,100)==0)
        fprintf(['Iter: ', num2str(i), ...
                 ' \t Avg. time/iter: ', num2str(tT/i*1000), ' ms \n']);
    end
end
tT = tT+toc(t0);
     
     
%% Errors                 
error_l = store_state(cord_nl([1,3,5]),:)'-phi_l_Traj;
error_r = store_state(cord_nl([2,4,6]),:)'-phi_r_Traj;
error_p = [vecnorm(error_l,2,2), vecnorm(error_r,2,2)];

eMAE = mean(abs([error_l error_r]));
eRMSE = sqrt(mean([error_l error_r].^2));

eRMSEp = mean([norm(eRMSE([1,3,5]),2) norm(eRMSE([2,4,6]),2)]);
eRMSEm = mean(eRMSE,2);
eMAEm  = mean(eMAE,2); % Old "avg_error"

fprintf(['-----------------------------------------\n', ...
         ' -- Total time: \t',num2str(tT),' s \n', ...
         ' -- Avg. t/iter:\t',num2str(tT/size(phi_l_Traj,1)*1000),' ms \n', ...
         ' -- Avg. error: \t',num2str(1000*eMAEm),' mm \n', ...
         ' -- Norm RMSE:  \t',num2str(1000*eRMSEp),' mm\n']);

     
%% PLOTS
time = 0:0.01:size(store_state,2)*0.01-0.01;

fig1 = figure(1);
fig1.Color = [1,1,1];
fig1.Units = 'normalized';
fig1.Position = [0.5 0.6 0.5 0.3];

subplot(7,2,1:2:12);
pa1som=plot(time,store_state(cord_nl([1 3 5]),:)','linewidth',1.5);
hold on
pa1ref=plot(time,phi_l_Traj(:,:),'--k','linewidth',1.2);
hold off
title('\textbf{Lower Left Corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 round(time(end))])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(7,2,2:2:12);
plot(time,store_state(cord_nl([2 4 6]),:)','linewidth',1.5)
hold on
plot(time,phi_r_Traj(:,:)','--k','linewidth',1.2)
hold off
title('\textbf{Lower Right Corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 round(time(end))])
set(gca, 'TickLabelInterpreter', 'latex');

Lgnd1 = legend([pa1som' pa1ref(1)], ...
               '$p_x\ (y_1,y_2)$','$p_y\ (y_3,y_4)$', '$p_z\ (y_5,y_6)$', '$r$', ...
               'Orientation','horizontal', 'Interpreter', 'latex');
Lgnd1.Position(1) = 0.5-Lgnd1.Position(3)/2;
Lgnd1.Position(2) = 0.02;
Lgnd1.Position(3:4) = Lgnd1.Position(3:4) + 0.01;

