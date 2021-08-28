% Calibration of the mass-spring-damper model taking as a reference the 
% inextensible cloth model of Franco et. al.
% Parameters to tune:
%   param.stiffness(1,3)
%   param.damping(1,3)
%   param.z_sum(1,1)
%
% Author: David Parent, davidparentalonso@gmail.com
% Last review: 10/02/2021

clear all; close all; clc

addpath('..\required_files\cloth_model_FColtraro')
addpath('..\required_files\cloth_model_DParent')
addpath('..\required_files\casadi-toolbox')
import casadi.*

opti = casadi.Opti();

% Define COM parameters(computation oriented model)
COM = struct;
COM.row = 4; COM.col = 4;
COM.mass = 0.1; 
COM.dt = 0.01;
COM.grav = 9.8;

% Define the SOM
nx = 10; ny = 10;
[SOM,pos] = initialize_model_parameters(nx, ny);

% Definite initial position of the nodes (defined here because we need it
% to compute ext_force
x_ini_SOM = [reshape(pos,[3*nx*ny 1]);zeros(3*nx*ny,1)]; %initial velocity=0
[reduced_pos,reduced_vel] = take_reduced_mesh(x_ini_SOM(1:3*nx*ny),x_ini_SOM((3*nx*ny+1):6*nx*ny));
x_ini_COM = [reduced_pos;reduced_vel];

% Define optimization variables
COM.z_sum = opti.variable(1,1);
COM.stiffness = opti.variable(1,3);
COM.damping = opti.variable(1,3);

COM.nodeInitial = lift_z(reshape(x_ini_COM(1:48), [16,3]), COM);

% Find initial spring length in each direction x,y,z
[mat_x, mat_y, mat_z] = compute_l0_linear(COM,1);
COM.mat_x = mat_x;
COM.mat_y = mat_y;
COM.mat_z = mat_z;

u_ini = x_ini_SOM(SOM.coord_controlados);
u_bef = u_ini;

nr = COM.row;
nc = COM.col;
cord_l = [1 nc 1+nr*nc nr*nc+nc 2*nr*nc+1 2*nr*nc+nc]; %coordenates of the linear model (of the lower corners)
cord_nl = [1 nx 1+nx*ny nx*ny+nx 2*nx*ny+1 2*nx*ny+nx]; %coordenates of the non-linear model (of the lower corners)

objectiu = 0;
for t=3:900
    % Enter the same "u_ini" to both models
    xnext_COM = linear_model_for_optimization(x_ini_COM,u_ini,u_bef,COM);
    xnext_SOM = cloth_simulator_firstorder(x_ini_SOM,u_ini,SOM);
    
    u_bef = u_ini;%update
    
    x_ini_COM = xnext_COM;
    x_ini_SOM = xnext_SOM(1:6*SOM.n_nodos);
      
    % Define a 3D trajectory to capture all behaviours
    if (t)>=0.0 && (t) <=50
        u_ini = u_ini;
    elseif (t)>50 && (t) <= 150
        u_ini(3) = u_ini(3) - 0.0005;
        u_ini(4) = u_ini(4) - 0.0005;
        u_ini(5) = u_ini(5) + 0.0005;
        u_ini(6) = u_ini(6) + 0.0005;
    elseif (t)>150 && (t) <= 350
        u_ini(1) = u_ini(1) - 0.0005;
        u_ini(2) = u_ini(2) - 0.0005;
        u_ini(3) = u_ini(3) - 0.001;
        u_ini(4) = u_ini(4) - 0.001;
    elseif (t)>350 && (t) <= 450
        u_ini(3) = u_ini(3) - 0.0005;
        u_ini(4) = u_ini(4) - 0.0005;
        u_ini(5) = u_ini(5) - 0.0005;
        u_ini(6) = u_ini(6) - 0.0005;
    elseif (t)>450 && (t) <= 550
        u_ini(1) = u_ini(1) + 0.0005;
        u_ini(2) = u_ini(2) + 0.0005;
        u_ini(3) = u_ini(3) + 0.0005;
        u_ini(4) = u_ini(4) + 0.0005;
    elseif (t) > 550 && (t) <= 750
        u_ini(1) = u_ini(1) + 0.001;
        u_ini(2) = u_ini(2) + 0.001;
        u_ini(3) = u_ini(3) + 0.001;
        u_ini(4) = u_ini(4) + 0.001;
        u_ini(5) = u_ini(5) - 0.0005;
        u_ini(6) = u_ini(6) - 0.0005;
    elseif (t)>750 && (t) <= 850
        u_ini(3) = u_ini(3) + 0.0005;
        u_ini(4) = u_ini(4) + 0.0005;
        u_ini(5) = u_ini(5) - 0.001;
        u_ini(6) = u_ini(6) - 0.001;
    else
        u_ini = u_ini;
    end
    
    objectiu = objectiu + (x_ini_COM(cord_l)-x_ini_SOM(cord_nl))'*(x_ini_COM(cord_l)-x_ini_SOM(cord_nl));
end

opti.minimize(objectiu);

p_opts = struct();
s_opts = struct('tol',10^-6,'max_iter',100);
opti.solver('ipopt',p_opts,s_opts);

% Constrains
opti.subject_to( COM.z_sum>0 ); 

% Initial seeds
opti.set_initial(COM.stiffness, [-300  -10 -200]);
opti.set_initial(COM.damping, [-4  -2   -4]);
opti.set_initial(COM.z_sum, 0.03);

sol = opti.solve();

% Print solution
sol.value(COM.stiffness)
sol.value(COM.damping)
sol.value(COM.z_sum)

% If optimization does not terminate, execute the following lines:
% opti.debug.value(COM.stiffness)
% opti.debug.value(COM.damping)
% opti.debug.value(COM.z_sum)
