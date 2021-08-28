% Solve the MPC problem where the COM is a linear cloth model and the SOM
% is the Coltraro's cloth model.
% Base code by: David Parent, davidparentalonso@gmail.com
% Modifications: Adrià Luque, adria.alados@gmail.com
% Last review: 24/08/2021

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
[A, B, f_ext] = create_model_linear_matrices(COM);

%% Start casADi optimization problem

% Lower corner coordinates
nr = COM.row; nc = COM.col;
coord_l =  [1 nc 1+nr*nc nr*nc+nc 2*nr*nc+1 2*nr*nc+nc];
coord_nl = [1 nx 1+nx*ny nx*ny+nx 2*nx*ny+1 2*nx*ny+nx];
COM.coord_controlados = [13 16 29 32 45 48];
COM.coord_lc = coord_l;


N = 30; % Prediction Horizon
ubound = 50e-3;

% Declare model variables
x = [SX.sym('pos',3*nr*nc,N+1);
     SX.sym('vel',3*nr*nc,N+1)];
u = SX.sym('u',6,N);
n_states = size(x,1); % 3*2*nxC*nyC

% Initial parameters of the optimization problem
P  = SX.sym('P', 1+6, max(n_states, N+1)); 
x0 = P(1, :)';
Rp = P(1+(1:6), 1:N+1);

% Optimization variables
w = u(:);
lbw = -ubound*ones(6*N,1);
ubw = +ubound*ones(6*N,1);

% Initialize other variables
objfun = 0;     % Objective function
g = [];         % Constraints
lbg = [];       % Lower bounds of g
ubg = [];       % Upper bounds of g

% Adaptive weight calculation:
%   From the current LCpos to the desired at the horizon
lc_dist = Rp(:,end) - x0(COM.coord_lc);
lc_dist = abs(lc_dist)/(norm(lc_dist)+eps);
Q = diag(lc_dist);
R = 20;

x(:,1) = x0;
for k = 1:N
    % Model Dynamics Constraint -> Definition
    x(:,k+1) = A*x(:,k) + B*u(:,k) + COM.dt*f_ext;

    % Objective function
    x_err = x(COM.coord_lc,k+1) - Rp(:,k+1);
    objfun = objfun + x_err'*Q*x_err + u(:,k)'*R*u(:,k);
end

opt_prob = struct('f', objfun, 'x', w, 'g', g, 'p', P);
config = struct;
config.print_time = 0;
config.ipopt.print_level = 0;  %0 min print - 3 max
config.ipopt.warm_start_init_point = 'yes'; %warm start

solver = nlpsol('ctrl_solver', 'ipopt', opt_prob, config);

%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP

%%
% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------

% Define trajectory to follow
Ref_l = load('required_files\trajectories\phi_l_3D.csv');
Ref_r = load('required_files\trajectories\phi_r_3D.csv');

u_ini = x_ini_SOM(SOM.coord_controlados);
u_bef = u_ini;

% Counters for the loop
sHp = 1; % start prediction horizon
Hp = N; % end prediction horizon

% Initialize things
in_params = zeros(1+6, max(n_states, N+1));
store_state(:,1) = x_ini_SOM;
store_state(:,2) = x_ini_SOM;
store_u(:,1) = zeros(6,1);
store_u(:,2) = zeros(6,1);

tT = 0;
t1 = tic;
%t0 = tic;
for i=3:size(Ref_l,1)
    Hp = Hp+1;
    sHp = sHp+1;

    % The last N+1 timesteps, trajectory should remain constant
    if i>=size(Ref_l,1)-(N+1)
        Ref_l_Hp = repmat(Ref_l(end,:), N+1,1);
        Ref_r_Hp = repmat(Ref_r(end,:), N+1,1);
    else
        Ref_l_Hp = Ref_l(sHp:Hp+1,:);
        Ref_r_Hp = Ref_r(sHp:Hp+1,:);
    end
    
    % Define input parameters for the optimizer (sliding window)
    in_params(1,:) = x_ini_COM';
    in_params(1+[1,3,5],1:N+1) = Ref_l_Hp';
    in_params(1+[2,4,6],1:N+1) = Ref_r_Hp';
    
    % Initial guess for optimizer (u: increments, guess UC=LC=Ref)
    args_x0 = [reshape(diff(in_params(1+(1:6),1:N),1,2),6*(N-1),1); zeros(6,1)];
    
    % Find the solution "sol"
    t0 = tic;
    sol = solver('x0', args_x0, 'lbx', lbw, 'ubx', ubw, ...
                 'lbg', lbg, 'ubg', ubg, 'p', in_params);
    
    u = reshape(full(sol.x),6,N)';
    u_lin = u(1,:)';
    u_SOM = u_lin + u_bef;
    u_bef = u_SOM;
    
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
tT1 = toc(t1);
     
%% Errors                 
error_l = store_state(coord_nl([1,3,5]),:)'-Ref_l;
error_r = store_state(coord_nl([2,4,6]),:)'-Ref_r;
error_p = [vecnorm(error_l,2,2), vecnorm(error_r,2,2)];

eMAE = mean(abs([error_l error_r]));
eRMSE = sqrt(mean([error_l error_r].^2));

eRMSEp = mean([norm(eRMSE([1,3,5]),2) norm(eRMSE([2,4,6]),2)]);
eRMSEm = mean(eRMSE,2);
eMAEm  = mean(eMAE,2); % Old "avg_error"

fprintf(['-----------------------------------------\n', ...
         ' -- Total time: \t',num2str(tT),' s \n', ...
         ' -- Avg. t/iter:\t',num2str(tT/size(Ref_l,1)*1000),' ms \n', ...
         ' -- Avg. error: \t',num2str(1000*eMAEm),' mm \n', ...
         ' -- Norm RMSE:  \t',num2str(1000*eRMSEp),' mm\n']);

     
%% PLOTS
time = 0:0.01:size(store_state,2)*0.01-0.01;

fig1 = figure(1);
fig1.Color = [1,1,1];
fig1.Units = 'normalized';
fig1.Position = [0.5 0.6 0.5 0.3];

subplot(7,2,1:2:12);
pa1som=plot(time,store_state(coord_nl([1 3 5]),:)','linewidth',1.5);
hold on
pa1ref=plot(time,Ref_l(:,:),'--k','linewidth',1.2);
hold off
title('\textbf{Lower Left Corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 round(time(end))])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(7,2,2:2:12);
plot(time,store_state(coord_nl([2 4 6]),:)','linewidth',1.5)
hold on
plot(time,Ref_r(:,:)','--k','linewidth',1.2)
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

