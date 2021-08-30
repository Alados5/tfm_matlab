function [Rwd, AllSt] = sim_ol_theta_realcloth(theta, opts)

addpath('..\required_files\cloth_model_New_L')

% Simulation Parameters
if nargin < 2
    % General Parameters
    nCOM = 4;
    NExp = 8;
    NTrial = 1;
    Ts = 0.01;
    ExpDate = '21_05_27';
    ExpSetN = 4;
else
    nCOM = opts.nCOM;
    NExp = opts.NExp;
    NTrial = opts.NTrial;
    Ts = opts.Ts;
    ExpDate = opts.ExpDate;
    ExpSetN = opts.ExpSetN;
end
% -------------------

% Define the SOM (REAL DATA)
SOMrec = load(['Gathered Data\Session_',ExpDate, ...
         '\Processed MAT S',num2str(ExpSetN),'\SOMtraj', ...
         num2str(NExp),'_',num2str(NTrial),'_',num2str(Ts*1000),'ms.mat']);
SOMrec = SOMrec.SOMrec;
SOMstate = SOMrec.states;
SOMtime = SOMrec.time;
nSOM = sqrt(size(SOMstate,1)/6);
Ts_rec = SOMtime(2) - SOMtime(1);
if (Ts-Ts_rec>eps)
    error("Input sample time (Ts) does not match the one on the recorded data"); 
end

SOM_node_ctrl = [nSOM^2-nSOM+1, nSOM^2];
SOM_coord_ctrl = [SOM_node_ctrl, ...
                  SOM_node_ctrl+nSOM^2, ...
                  SOM_node_ctrl+2*nSOM^2];

% Input trajectory are the real upper corners
TrajU = SOMstate(SOM_coord_ctrl,:);
u_ini = TrajU(:,1);


% Define COM parameters
COM = struct;
COM.row = nCOM;
COM.col = nCOM;
COM.mass = 0.1;
COM.grav = 9.8;
COM.dt = Ts;

% Parameters: theta
COM.stiffness = theta.stiffness; 
COM.damping = theta.damping; 
COM.z_sum = theta.z_sum;

% Controlled coordinates (upper corners in x,y,z)
COM_node_ctrl = [nCOM^2-nCOM+1, nCOM^2];
COM.coord_ctrl = [COM_node_ctrl, ...
                  COM_node_ctrl+nCOM^2, ...
                  COM_node_ctrl+2*nCOM^2];

% Define initial state of the nodes (needed for ext_force)
x_ini_SOM = SOMstate(:, 1);
[pos_rd,~] = take_reduced_mesh(x_ini_SOM(1:3*nSOM^2),x_ini_SOM(3*nSOM^2+1:6*nSOM^2), nSOM, nCOM);
x_ini_COM = [pos_rd; zeros(3*nCOM^2,1)];

% Initial position of the nodes
COM.nodeInitial = lift_z(reshape(x_ini_COM(1:3*nCOM^2), [nCOM^2,3]), COM);

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
    
    pos_SOMi = SOMstate(1:3*nSOM^2, tk);
    vel_SOMi = SOMstate(3*nSOM^2+1:6*nSOM^2, tk);
    SOMstv(:,tk) = [pos_SOMi; vel_SOMi];
    
    [pos_rd, vel_rd] = take_reduced_mesh(pos_SOMi,vel_SOMi, nSOM, nCOM);
    SOMRstv(:,tk) = [pos_rd; vel_rd];
    
    u_bef = u_SOM; % Update

end

% Convert to cm and square to penalize big differences more
avg_lin_error = mean((100*(SOMRstv-COMstv)).^2,2);
avg_lin_error_pos = avg_lin_error(1:3*nCOM^2);

% Ponderate to penalize lower corners more
err_mask = kron([1 1 1]', (floor(nCOM-1/nCOM:-1/nCOM:0)'+1)/nCOM);
wavg_lin_error_pos = avg_lin_error_pos.*err_mask.^2;

% Final Reward
Rwd = -norm(wavg_lin_error_pos, 1);

fprintf([' Reward: ', num2str(Rwd), '\n']);

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
    scatter3(AllSt.SOM(1:nSOM^2,tk), ...
             AllSt.SOM(nSOM^2+1:2*nSOM^2,tk), ...
             AllSt.SOM(2*nSOM^2+1:3*nSOM^2,tk));
    hold on
    scatter3(AllSt.COM(1:nCOM^2,tk), ...
             AllSt.COM(nCOM^2+1:2*nCOM^2,tk), ...
             AllSt.COM(2*nCOM^2+1:3*nCOM^2,tk));
    hold off
    axis equal
    pause(1e-6);
end
%}
%{
coord_ctrl = COM.coord_ctrl;
coord_lc = [1 nCOM 1+nCOM^2 nCOM+nCOM^2 1+2*nCOM^2 nCOM+2*nCOM^2];
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
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
%}

end