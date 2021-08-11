% Calibration of the mass-spring-damper model taking as a reference the 
% inextensible cloth model of Franco et. al. (old version)
% Parameters to tune:
%   param.stiffness(1,3)
%   param.damping(1,3)
%   param.z_sum(1,1)

clear; close all; clc

addpath('..\required_files\cloth_model_FColtraro')
addpath('..\required_files\cloth_model_DParent')
addpath('..\required_files\casadi-toolbox')
import casadi.*

opti = casadi.Opti();

SOMsz = 10;
COMsz = 4;
Ts = 0.005;

% Define COM parameters
COM = struct;
COM.row = COMsz;
COM.col = COMsz;
COM.mass = 0.1; 
COM.dt = Ts;
COM.grav = 9.8;

% Block sizes
COMlength = COM.row*COM.col;
blksz = 3*COMlength;

% Controlled coordinates: upper corners
COM_nodes_ctrl = [COMsz*(COMsz-1)+1, COMsz*COMsz];
COM.coord_ctrl = [COM_nodes_ctrl, COM_nodes_ctrl+COMsz^2, COM_nodes_ctrl+2*COMsz^2];

% Define the SOM
nx = SOMsz;
ny = SOMsz;
[SOM,pos] = initialize_model_parameters(nx, ny, Ts);

% Definite initial position of the nodes (defined here because we need it
% to compute ext_force
x_ini_SOM = [reshape(pos,[3*nx*ny 1]);zeros(3*nx*ny,1)]; %initial velocity=0
[reduced_pos,reduced_vel] = take_reduced_mesh(x_ini_SOM(1:3*nx*ny), ...
                                              x_ini_SOM((3*nx*ny+1):6*nx*ny), ...
                                              SOMsz, COMsz);
x_ini_COM = [reduced_pos; reduced_vel];

% Define optimization variables
COM.z_sum = opti.variable(1,1);
COM.stiffness = opti.variable(1,3);
COM.damping = opti.variable(1,3);

COM.nodeInitial = lift_z(reshape(x_ini_COM(1:blksz), [COMlength,3]), COM);

% Find initial spring length in each direction x,y,z
[mat_x, mat_y, mat_z] = compute_l0_linear(COM,1);
COM.mat_x = mat_x;
COM.mat_y = mat_y;
COM.mat_z = mat_z;

u_ini = x_ini_SOM(SOM.coord_controlados);
u_bef = u_ini;
u_SOM = u_ini;

nr = COM.row;
nc = COM.col;

% Lower corner coordinates for both linear and nonlinear models
coord_l  = [1 nc 1+nr*nc nr*nc+nc 2*nr*nc+1 2*nr*nc+nc];
coord_nl = [1 nx 1+nx*ny nx*ny+nx 2*nx*ny+1 2*nx*ny+nx];

NT = 50;

objectiu = 0;
for t=3:18*NT
    % Enter the same "u_ini" to both models
    u_lin = u_SOM - u_bef;
    xnext_COM = linear_model_for_optimization(x_ini_COM,u_lin,COM);
    xnext_SOM = cloth_simulator_firstorder(x_ini_SOM,u_SOM,SOM);
    
    u_bef = u_SOM; % Update
    
    x_ini_COM = xnext_COM;
    x_ini_SOM = xnext_SOM(1:6*SOM.n_nodos);
    
    [reduced_pos,reduced_vel] = take_reduced_mesh(x_ini_SOM(1:3*nx*ny), ...
                                                  x_ini_SOM((3*nx*ny+1):6*nx*ny), ...
                                                  SOMsz, COMsz);
    x_red_SOM = [reduced_pos; reduced_vel];
      
    % Define a 3D trajectory to capture all behaviours
    if t <= NT
        u_SOM = u_ini;
    elseif t <= 3*NT
        u_SOM(3) = u_SOM(3) - 0.0005;
        u_SOM(4) = u_SOM(4) - 0.0005;
        u_SOM(5) = u_SOM(5) + 0.0005;
        u_SOM(6) = u_SOM(6) + 0.0005;
    elseif t <= 7*NT
        u_SOM(1) = u_SOM(1) - 0.0005;
        u_SOM(2) = u_SOM(2) - 0.0005;
        u_SOM(3) = u_SOM(3) - 0.001;
        u_SOM(4) = u_SOM(4) - 0.001;
    elseif t <= 9*NT
        u_SOM(3) = u_SOM(3) - 0.0005;
        u_SOM(4) = u_SOM(4) - 0.0005;
        u_SOM(5) = u_SOM(5) - 0.0005;
        u_SOM(6) = u_SOM(6) - 0.0005;
    elseif t <= 11*NT
        u_SOM(1) = u_SOM(1) + 0.0005;
        u_SOM(2) = u_SOM(2) + 0.0005;
        u_SOM(3) = u_SOM(3) + 0.0005;
        u_SOM(4) = u_SOM(4) + 0.0005;
    elseif t <= 15*NT
        u_SOM(1) = u_SOM(1) + 0.001;
        u_SOM(2) = u_SOM(2) + 0.001;
        u_SOM(3) = u_SOM(3) + 0.001;
        u_SOM(4) = u_SOM(4) + 0.001;
        u_SOM(5) = u_SOM(5) - 0.0005;
        u_SOM(6) = u_SOM(6) - 0.0005;
    %{x
    elseif t <= 17*NT
        u_SOM(3) = u_SOM(3) + 0.0005;
        u_SOM(4) = u_SOM(4) + 0.0005;
        u_SOM(5) = u_SOM(5) - 0.001;
        u_SOM(6) = u_SOM(6) - 0.001;
    %}
    else
        u_SOM = u_SOM;
    end
    
    objectiu = objectiu + (x_ini_COM(coord_l)-x_ini_SOM(coord_nl))'*(x_ini_COM(coord_l)-x_ini_SOM(coord_nl));
    %objectiu = objectiu + (x_ini_COM-x_red_SOM)' * (x_ini_COM-x_red_SOM);
end

opti.minimize(objectiu);

% Maximum Iterations
MaxIter = 200;

p_opts = struct();
s_opts = struct('tol',1e-6, 'max_iter',MaxIter);
opti.solver('ipopt', p_opts, s_opts);

% Constrains
opti.subject_to( COM.z_sum>=0 ); 

% Initial seeds
%{-
opti.set_initial(COM.stiffness, [-300 -10 -200]);
opti.set_initial(COM.damping, [-4 -2 -4]);
opti.set_initial(COM.z_sum, 0.03);
%}
%{
opti.set_initial(COM.stiffness, [-300 -10 -400]);
opti.set_initial(COM.damping, [-4 -2 -5]);
opti.set_initial(COM.z_sum, 0.04);
%}
sol = opti.solve();

% Print solution
sol.value(COM.stiffness)
sol.value(COM.damping)
sol.value(COM.z_sum)

% If optimization does not terminate, execute the following lines:
%{
opti.debug.value(COM.stiffness)
opti.debug.value(COM.damping)
opti.debug.value(COM.z_sum)
%}



