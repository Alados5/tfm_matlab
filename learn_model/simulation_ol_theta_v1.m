function [Rwd, AllSt] = simulation_ol_theta_v1(theta, opts)

addpath('..\required_files\cloth_model_FColtraro')
addpath('..\required_files\cloth_model_DParent')

% Simulation Parameters
if nargin < 2
    % General Parameters
    nCOM = 4;
    nSOM = 10;
    Ts = 0.010;
    NExp = 4;
else
    nCOM = opts.nCOM;
    nSOM = opts.nSOM;
    Ts = opts.Ts; 
    NExp = opts.NExp;
end
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

% Parameters: theta
COM.stiffness = theta.stiffness; 
COM.damping = theta.damping; 
COM.z_sum = theta.z_sum;


% Controlled coordinates (upper corners in x,y,z)
COM_node_ctrl = [nxC*(nyC-1)+1, nxC*nyC];
COM.coord_ctrl = [COM_node_ctrl, ...
                  COM_node_ctrl+nxC*nyC, ...
                  COM_node_ctrl+2*nxC*nyC];

% Define the SOM (NONLINEAR)
nxS = nSOM;
nyS = nSOM;
SOMlength = nxS*nyS;
[SOM, pos] = initialize_model_parameters(nxS, nyS, Ts);

% Controlled coordinates (upper corners in x,y,z)
SOM_node_ctrl = [nxS*(nyS-1)+1, nxS*nyS];
SOM.coord_ctrl = [SOM_node_ctrl SOM_node_ctrl+nxS*nyS SOM_node_ctrl+2*nxS*nyS];

% Define initial position of the nodes (needed for ext_force)
% Second half is velocity (initial v=0)
x_ini_SOM = [reshape(pos,[3*nxS*nyS 1]); zeros(3*nxS*nyS,1)];
[pos_rd,~] = take_reduced_mesh(x_ini_SOM(1:3*nxS*nyS),x_ini_SOM(3*nxS*nyS+1:6*nxS*nyS), nSOM, nCOM);
x_ini_COM = [pos_rd; zeros(3*nxC*nyC,1)];

% Initial position of the nodes
COM.nodeInitial = lift_z(reshape(x_ini_COM(1:3*nxC*nyC), [nxC*nyC,3]), COM);


% Find initial spring length in each direction x,y,z
[COM.mat_x, COM.mat_y, COM.mat_z] = compute_l0_linear(COM,0);

% Find linear matrices
[A_COM, B_COM, f_COM] = create_model_linear_matrices(COM);


% Define a rich trajectory for the upper corners
TrajU = load(['TrajUs/TrajU_',num2str(NExp),'.mat']);
TrajU = TrajU.TrajU;

% Either initial state or input (same for v1)
u_ini = x_ini_SOM(SOM.coord_controlados);
u_bef = u_ini;

SOMstv(:,1:2) = [x_ini_SOM x_ini_SOM];
COMstv(:,1:2) = [x_ini_COM x_ini_COM];
SOMRstv(:,1:2) = [x_ini_COM x_ini_COM];
for tk=3:size(TrajU,2)
    
    u_SOM = TrajU(:,tk);
    u_lin = u_SOM - u_bef;
    
    xnext_COM = A_COM*COMstv(:,tk-1) + B_COM*u_lin + Ts*f_COM;
    COMstv(:,tk) = xnext_COM;
    
    xnext_SOM = cloth_simulator_secondorder([SOMstv(:,tk-1);SOMstv(:,tk-2)],u_SOM,SOM);
    SOMstv(:,tk) = xnext_SOM(1:6*SOM.n_nodos);
    
    pos_SOMi = SOMstv(1:3*SOMlength, tk);
    vel_SOMi = SOMstv((1+3*SOMlength):6*SOMlength, tk);
    [pos_rd, vel_rd] = take_reduced_mesh(pos_SOMi,vel_SOMi, nSOM, nCOM);
    
    x_red_SOM = [pos_rd; vel_rd];
    SOMRstv(:,tk) = x_red_SOM;
    
    u_bef = u_SOM; % Update

end


avg_lin_error = mean(abs(SOMRstv-COMstv),2);
avg_lin_error_pos = avg_lin_error(1:3*COMlength);
Rwd = -norm(avg_lin_error_pos, 1);

disp(' ');
disp(['Reward: ', num2str(Rwd)]);

if isnan(Rwd)
    Rwd = -inf;
end

AllSt = struct();
AllSt.SOM = SOMstv;
AllSt.SOMrd = SOMRstv;
AllSt.COM = COMstv;
AllSt.uSOM = TrajU;

end