function [Rwd, AllSt] = sim_ol_theta_nlmdl_new(theta, opts)

addpath('..\required_files\cloth_model_New_NL')
addpath('..\required_files\cloth_model_New_L')

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


% Load a rich trajectory for the upper corners
TrajU = load(['TrajUs/TrajU_',num2str(NExp),'.mat']);
TrajU = TrajU.TrajU;
u_ini = TrajU(:,1);

lCloth = norm(u_ini([2 4 6]) - u_ini([1 3 5]));
cCloth = (u_ini([2 4 6]) + u_ini([1 3 5]))/2 - [0; 0; lCloth/2];


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
%SOMlength = nxS*nyS;
[SOM, pos] = initialize_nl_model(lCloth,nSOM,cCloth,Ts);

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


% Initialize variables
u_bef = u_ini;
SOMstv(:,1) = x_ini_SOM;
COMstv(:,1) = x_ini_COM;
SOMRstv(:,1) = x_ini_COM;

% Iterate (Open-loop)
for tk=2:size(TrajU,2)
    
    u_SOM = TrajU(:,tk);
    u_lin = u_SOM - u_bef;
    
    x_COMi = A_COM*COMstv(:,tk-1) + B_COM*u_lin + Ts*f_COM;
    COMstv(:,tk) = x_COMi;
    
    [pos_SOMi, vel_SOMi] = simulate_cloth_step(SOMstv(:,tk-1), u_SOM, SOM);
    SOMstv(:,tk) = [pos_SOMi; vel_SOMi];
    
    [pos_rd, vel_rd] = take_reduced_mesh(pos_SOMi,vel_SOMi, nSOM, nCOM);
    SOMRstv(:,tk) = [pos_rd; vel_rd];
    
    u_bef = u_SOM; % Update

end

% Convert to cm and square to penalize big differences more
avg_lin_error = mean((100*(SOMRstv-COMstv)).^2,2);
avg_lin_error_pos = avg_lin_error(1:3*COMlength);

% Ponderate to penalize lower corners more
err_mask = kron([1 1 1]', (floor(nCOM-1/nCOM:-1/nCOM:0)'+1)/nCOM);
wavg_lin_error_pos = avg_lin_error_pos.*err_mask.^2;

% Final Reward
Rwd = -norm(wavg_lin_error_pos, 1);

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

%{
for tk=1:size(TrajU,2)
    scatter3(AllSt.SOM(1:nxS*nyS,tk), ...
             AllSt.SOM(nxS*nyS+1:2*nxS*nyS,tk), ...
             AllSt.SOM(2*nxS*nyS+1:3*nxS*nyS,tk));
    axis equal
    pause(1e-6);
end
%}
%{
plot3(AllSt.SOMrd(coord_ctrl(1),:), AllSt.SOMrd(coord_ctrl(3),:), AllSt.SOMrd(coord_ctrl(5),:), '-k')
hold on
plot3(AllSt.SOMrd(coord_ctrl(2),:), AllSt.SOMrd(coord_ctrl(4),:), AllSt.SOMrd(coord_ctrl(6),:), '-k')
plot3(AllSt.SOMrd(coord_lc(1),:), AllSt.SOMrd(coord_lc(3),:), AllSt.SOMrd(coord_lc(5),:), '-b')
plot3(AllSt.SOMrd(coord_lc(2),:), AllSt.SOMrd(coord_lc(4),:), AllSt.SOMrd(coord_lc(6),:), '-b')
plot3(AllSt.COM(coord_lc(1),:), AllSt.COM(coord_lc(3),:), AllSt.COM(coord_lc(5),:), '--r')
plot3(AllSt.COM(coord_lc(2),:), AllSt.COM(coord_lc(4),:), AllSt.COM(coord_lc(6),:), '--r')
hold off
box on
grid on
axis equal
%}

end