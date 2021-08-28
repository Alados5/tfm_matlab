clear; close all; clc;

%% Parameters

addpath('required_files\cloth_model_New_L')
addpath('required_files\cloth_model_New_NL')
addpath('..\..\required_files\casadi-toolbox')
import casadi.*

% Unstabilized with Ts=0.012 and step_t=1 (y-movement)

% General Parameters
Ts = 0.010;
nSOM = 10;
nCOM = 4;

lCloth = 0.3;
cCloth = [0 -0.4 0.1];
aCloth = 0;


%% Input step defnition

%step_p = 0.1*[1; 1; 0; 0; 0; 0]; %m
%step_p = 0.1*[0; 0; -1; -1; 0; 0]; %m
%step_p = 0.1*[0; 0; 0; 0; 1; 1]; %m
%step_p =  0.1*[1; -1; 0; 0; 0; 0]; %m
step_p = 0.1*[1; 1; -3; -3; 2; 2]; %m
%{
% Generic case
step_t = 0.1; %s
u_traj = [zeros(6,round(0.5/Ts)), ...
          step_p/step_t*Ts.*ones(6,round(step_t/Ts)), ...
          zeros(6,round(2.5/Ts))];
%}
%{x
% Progressive acceleration
step_t = 0.4; %s
u_traj = [zeros(6,round(0.5/Ts)), ...
          step_p/step_t*Ts.*ones(6,round(step_t/Ts)), ...
        2*step_p/step_t*Ts.*ones(6,round(step_t/Ts)), ...
        4*step_p/step_t*Ts.*ones(6,round(step_t/Ts)), ...
        2*step_p/step_t*Ts.*ones(6,round(step_t/Ts)), ...
          step_p/step_t*Ts.*ones(6,round(step_t/Ts)), ...
          zeros(6,round(1/Ts))];
%}
%{
% For size 7x7
step_p = 0.1*[0; 0; -1; -1; 0; 0]; %m
step_t = 10*0.1; %s
u_traj = [zeros(6,round(0.5/Ts)), ...
          step_p/step_t*Ts.*ones(6,round(step_t/Ts)), ...
          zeros(6,round(0.5/Ts))];
%}
nPtTraj = size(u_traj,2);


%% Model Initialization

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
COM.stiffness = [-305.6028  -13.4221 -225.8987]; 
COM.damping = [-4.0042   -2.5735   -3.9090]; 
COM.z_sum = 0.0312;

% Important Coordinates (upper and lower corners in x,y,z)
COM_nd_ctrl = [nxC*(nyC-1)+1, nxC*nyC];
COM.coord_ctrl = [COM_nd_ctrl, COM_nd_ctrl+nxC*nyC, COM_nd_ctrl+2*nxC*nyC];
C_coord_lc = [1 nyC 1+nxC*nyC nxC*nyC+nyC 2*nxC*nyC+1 2*nxC*nyC+nyC]; 
COM.coord_lc = C_coord_lc;


% Define the SOM (NONLINEAR)
nxS = nSOM;
nyS = nSOM;
SOMlength = nxS*nyS;
[SOM, pos] = initialize_nl_model(lCloth,nSOM,cCloth,aCloth,Ts);
S_coord_lc = SOM.coord_lc;


% Define initial position of the nodes (needed for ext_force)
% Second half is velocity (initial v=0)
x_ini_SOM = [reshape(pos,[3*nxS*nyS 1]); zeros(3*nxS*nyS,1)];

% Reduce initial SOM position to COM size if necessary
[reduced_pos,~] = take_reduced_mesh(x_ini_SOM(1:3*nxS*nyS),x_ini_SOM(3*nxS*nyS+1:6*nxS*nyS), nSOM, nCOM);
x_ini_COM = [reduced_pos; zeros(3*nxC*nyC,1)];

% Rotate initial COM position to XZ plane
RCloth_ini = [cos(aCloth) -sin(aCloth) 0; sin(aCloth) cos(aCloth) 0; 0 0 1];
posCOM = reshape(x_ini_COM(1:3*nxC*nyC), [nxC*nyC,3]);
posCOM_XZ = (RCloth_ini^-1 * posCOM')';

% Initial position of the nodes
COM.nodeInitial = lift_z(posCOM_XZ, COM);

% Find initial spring length in each direction x,y,z
[COM.mat_x, COM.mat_y, COM.mat_z] = compute_l0_linear(COM,0);

% Find linear matrices
[A_COM, B_COM, f_COM] = create_model_linear_matrices(COM);


% Initialize control
u_ini  = x_ini_SOM(SOM.coord_ctrl);
u_bef  = u_ini;
u_SOM  = u_ini;

% Simulate some SOM steps to stabilize the NL model
warning('off','MATLAB:nearlySingularMatrix');
lastwarn('','');
[p_ini_SOM, ~] = simulate_cloth_step(x_ini_SOM,u_SOM,SOM);
[~, warnID] = lastwarn;
while strcmp(warnID, 'MATLAB:nearlySingularMatrix')
    lastwarn('','');
    x_ini_SOM = [p_ini_SOM; zeros(3*nxS*nyS,1)];
    [p_ini_SOM, ~] = simulate_cloth_step(x_ini_SOM,u_SOM,SOM);
    [~, warnID] = lastwarn;
end
warning('on','MATLAB:nearlySingularMatrix');


%% Simulation

x_COM = x_ini_COM;
x_SOM = x_ini_SOM;
p_SOM = x_SOM(1:3*SOMlength); v_SOM = x_SOM(3*SOMlength+1:end);
[p_SOMr, v_SOMr] = take_reduced_mesh(p_SOM, v_SOM, nSOM, nCOM);
x_SOMr = [p_SOMr; v_SOMr];

states_lin = zeros(6*nxC*nyC,nPtTraj+1);
states_nlr = zeros(6*nxC*nyC,nPtTraj+1);
states_lin(:,1) = x_COM;
states_nlr(:,1) = x_SOMr;
time = 0:Ts:Ts*nPtTraj;

for tk=1:nPtTraj
    
    u_lin = u_traj(:,tk);
    u_SOM = u_lin + u_bef;
    u_bef = u_SOM;
    
    x_COM = A_COM*x_COM + B_COM*u_lin + COM.dt*f_COM;
    [p_SOM, v_SOM] = simulate_cloth_step(x_SOM,u_SOM,SOM);
    x_SOM = [p_SOM; v_SOM];
    
    [p_SOMr, v_SOMr] = take_reduced_mesh(p_SOM, v_SOM, nSOM, nCOM);
    x_SOMr = [p_SOMr; v_SOMr];
    
    states_lin(:,tk+1) = x_COM;
    states_nlr(:,tk+1) = x_SOMr;
    
end


%% Plots

fig1 = figure(1);
fig1.Units = 'normalized';
fig1.Position = [0.4 0.4 0.6 0.4];
fig1.Color = [1,1,1];

subplot(10,2,9:2:18);
pa_nlr=plot(time, states_nlr(COM.coord_lc([1 3 5]),:),'linewidth',1.5);
hold on
pa_lin=plot(time, states_lin(COM.coord_lc([1 3 5]),:),'--','linewidth',1.5);
hold off
title('\textbf{Lower Left Corner Response}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex','FontSize',10);
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',10);
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex','FontSize',10);
ylim_l = ylim;

subplot(10,2,10:2:18);
plot(time, states_nlr(COM.coord_lc([2 4 6]),:),'linewidth',1.5)
hold on
plot(time, states_lin(COM.coord_lc([2 4 6]),:),'--','linewidth',1.5)
hold off
title('\textbf{Lower Right Corner Response}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex','FontSize',10);
ylabel('Position [m]', 'Interpreter', 'latex','FontSize',10);
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex','FontSize',10);
ylim_r = ylim;

ylim_c = [min(ylim_l(1), ylim_r(1)), max(ylim_l(2), ylim_r(2))];
subplot(10,2,9:2:18);
ylim(ylim_c);
subplot(10,2,10:2:18);
ylim(ylim_c);

subplot(10,2,1:2:4);
pa_u1 = plot(time(1:end-1), u_traj(1,:),'b','linewidth',1);
hold on
pa_u3 = plot(time(1:end-1), u_traj(3,:),'r','linewidth',1);
pa_u5 = plot(time(1:end-1), u_traj(5,:),'color',[1,0.6,0],'linewidth',1);
hold off
title('\textbf{Upper Left Corner Input}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex','FontSize',10);
ylabel('$u$ [m]', 'Interpreter', 'latex','FontSize',10);
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex','FontSize',10);

subplot(10,2,2:2:4);
plot(time(1:end-1), u_traj(2,:),'b','linewidth',1)
hold on
plot(time(1:end-1), u_traj(4,:),'r','linewidth',1)
plot(time(1:end-1), u_traj(6,:),'color',[1,0.6,0],'linewidth',1)
title('\textbf{Upper Right Corner Input}', 'Interpreter', 'latex')
hold off
grid on
xlabel('Time [s]', 'Interpreter', 'latex','FontSize',10);
ylabel('$u$ [m]', 'Interpreter', 'latex','FontSize',10);
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex','FontSize',10);


Lgnd1 = legend([pa_nlr', pa_lin', pa_u1, pa_u3, pa_u5], ...
               '$y_{1,2}^{NL}$','$y_{3,4}^{NL}$','$y_{5,6}^{NL}$ \quad', ...
               '$y_{1,2}^{L}$', '$y_{3,4}^{L}$', '$y_{5,6}^{L}$  \quad', ...
               '$u_{1,2}$',     '$u_{3,4}$',     '$u_{5,6}$',...
               'Orientation','horizontal', 'Interpreter', 'latex');
Lgnd1.Position(1) = 0.5-Lgnd1.Position(3)/2;
Lgnd1.Position(2) = 0.02;
Lgnd1.Position(3:4) = Lgnd1.Position(3:4) + 0.01;









