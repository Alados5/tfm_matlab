% Solve the MPC problem where the COM is a linear cloth model and the SOM
% is the Coltraro's cloth model.
% Base code by: David Parent, davidparentalonso@gmail.com
% Modifications: Adrià Luque, adria.alados@gmail.com
% Last review: 24/08/2021

clear; close all; clc;

Ts = 0.01;
NHp = 30; % Prediction Horizon
ubound = 50e-3;
nCOM = 4;
nSOM = 10;
Qk = 1;
R = 20;

rot_traj = 0;
rot_pars = 0;

addpath('required_files\cloth_model_FColtraro')
addpath('required_files\cloth_model_New_L')
addpath('..\..\required_files\casadi-toolbox')
import casadi.*

% Define trajectory to follow
RTraj = '3D';
Ref_l = load(['required_files\trajectories\phi_l_',RTraj,'.csv']);
Ref_r = load(['required_files\trajectories\phi_r_',RTraj,'.csv']);

% ROT90 CHANGE: Rotate ref 90 degrees
if rot_traj == 1
    Ref_aux = Ref_l;
    Ref_l(:,1) =  Ref_aux(:,2);
    Ref_l(:,2) = -Ref_aux(:,1);
    Ref_aux = Ref_r;
    Ref_r(:,1) =  Ref_aux(:,2);
    Ref_r(:,2) = -Ref_aux(:,1);
end

% Define COM parameters(computation oriented model)
COM = struct;
COM.row = nCOM;
COM.col = nCOM;
COM.mass = 0.1;
COM.grav = 9.8;
COM.dt = Ts;

% parameters from calibration
COM.stiffness = [-305.6028  -13.4221 -225.8987]; 
COM.damping = [-4.0042   -2.5735   -3.9090]; 
COM.z_sum = 0.0312;

% ROT90 CHANGE: Swap x/y parameters
if rot_pars == 1
    COM.stiffness(1:2) = COM.stiffness([2,1]);
    COM.damping(1:2) = COM.damping([2,1]);
end


% Define the SOM
nx = nSOM; %COLS!
ny = nSOM; %ROWS!
[SOM, pos] = initialize_model_parameters(nx,ny,Ts);

% ROT90 CHANGE: Rotate pos 90 degrees
if rot_traj == 1
    posaux = pos;
    pos(:,1) =  posaux(:,2);
    pos(:,2) = -posaux(:,1);
end


% Important Coordinates (upper and lower corners in x,y,z)
nr = COM.row;
nc = COM.col;
coord_l =  [1 nc 1+nr*nc nr*nc+nc 2*nr*nc+1 2*nr*nc+nc];
coord_nl = [1 nx 1+nx*ny nx*ny+nx 2*nx*ny+1 2*nx*ny+nx];
COM_nd_ctrl = [nr*(nc-1)+1, nr*nc];
COM.coord_ctrl = [COM_nd_ctrl, COM_nd_ctrl+nr*nc, COM_nd_ctrl+2*nr*nc];
COM.coord_lc = coord_l;

% Define initial position of the nodes (needed for ext_force)
% Second half is velocity (initial v=0)
Nd_C = nr*nc;
Nd_S = nx*ny;
x_ini_SOM = [reshape(pos,[3*Nd_S 1]);zeros(3*Nd_S,1)]; %v0=0
[reduced_pos,~] = take_reduced_mesh(x_ini_SOM(1:3*Nd_S),x_ini_SOM(3*Nd_S+1:6*Nd_S));
x_ini_COM = [reduced_pos;zeros(3*Nd_C,1)];

% Initial position of the COM nodes
COM.nodeInitial = lift_z(reshape(x_ini_COM(1:3*Nd_C), [Nd_C,3]), COM);

% Find initial spring length in each direction x,y,z
[COM.mat_x, COM.mat_y, COM.mat_z] = compute_l0_linear(COM,0);

% Find linear matrices
[A, B, f_ext] = create_model_linear_matrices(COM);

%% Start casADi optimization problem

% Declare model variables
x = [SX.sym('pos',3*nr*nc,NHp+1);
     SX.sym('vel',3*nr*nc,NHp+1)];
u = SX.sym('u',6,NHp);
n_states = size(x,1); % 3*2*nxC*nyC

% Initial parameters of the optimization problem
P  = SX.sym('P', 1+6, max(n_states, NHp+1)); 
x0 = P(1, :)';
Rp = P(1+(1:6), 1:NHp+1);

% Optimization variables
w = u(:);
lbw = -ubound*ones(6*NHp,1);
ubw = +ubound*ones(6*NHp,1);

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

x(:,1) = x0;
for k = 1:NHp
    % Model Dynamics Constraint -> Definition
    x(:,k+1) = A*x(:,k) + B*u(:,k) + COM.dt*f_ext;

    % Objective function
    x_err = x(COM.coord_lc,k+1) - Rp(:,k+1);
    objfun = objfun + x_err'*Qk*Q*x_err + u(:,k)'*R*u(:,k);
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

u_ini = x_ini_SOM(SOM.coord_controlados);
u_bef = u_ini;

% Counters for the loop
sHp = 1; % start prediction horizon
Hp = NHp; % end prediction horizon

% Initialize things
in_params = zeros(1+6, max(n_states, NHp+1));
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
    if i>=size(Ref_l,1)-(NHp+1)
        Ref_l_Hp = repmat(Ref_l(end,:), NHp+1,1);
        Ref_r_Hp = repmat(Ref_r(end,:), NHp+1,1);
    else
        Ref_l_Hp = Ref_l(sHp:Hp+1,:);
        Ref_r_Hp = Ref_r(sHp:Hp+1,:);
    end
    
    % Define input parameters for the optimizer (sliding window)
    in_params(1,:) = x_ini_COM';
    in_params(1+[1,3,5],1:NHp+1) = Ref_l_Hp';
    in_params(1+[2,4,6],1:NHp+1) = Ref_r_Hp';
    
    % Initial guess for optimizer (u: increments, guess UC=LC=Ref)
    args_x0 = [reshape(diff(in_params(1+(1:6),1:NHp),1,2),6*(NHp-1),1); zeros(6,1)];
    
    % Find the solution "sol"
    t0 = tic;
    sol = solver('x0', args_x0, 'lbx', lbw, 'ubx', ubw, ...
                 'lbg', lbg, 'ubg', ubg, 'p', in_params);
    
    u = reshape(full(sol.x),6,NHp)';
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
time = 0:Ts:size(store_state,2)*Ts-Ts;

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


%% 3D PLOT

fig3 = figure(3);
fig3.Color = [1,1,1];

SOM_lowc = [1 ny]; 
SOM_nd_ctrl = [nx*(ny-1)+1, nx*ny];
SOM.coord_ctrl = [SOM_nd_ctrl, SOM_nd_ctrl+nx*ny, SOM_nd_ctrl+2*nx*ny];

store_pos = store_state(1:3*Nd_S,:);
store_x = store_pos(1:Nd_S,:);
limx = [floor(min(store_x(:))*10), ceil(max(store_x(:))*10)]/10;
store_y = store_pos(Nd_S+1:2*Nd_S,:);
limy = [floor(min(store_y(:))*10), ceil(max(store_y(:))*10)]/10;
store_z = store_pos(2*Nd_S+1:3*Nd_S,:);
limz = [floor(min(store_z(:))*10), ceil(max(store_z(:))*10)]/10;

All_uSOM = store_state(SOM.coord_ctrl,:);

plot3(All_uSOM(1:2,:)',All_uSOM(3:4,:)',All_uSOM(5:6,:)');
hold on
plot3(store_x(SOM_lowc,:)', store_y(SOM_lowc,:)', store_z(SOM_lowc,:)');
hold off
axis equal; box on; grid on;
xlim(limx);
ylim(limy);
zlim(limz);
set(gca, 'TickLabelInterpreter','latex');
xlabel('X', 'Interpreter','latex');
ylabel('Y', 'Interpreter','latex');
zlabel('Z', 'Interpreter','latex');




