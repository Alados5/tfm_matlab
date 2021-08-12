%{
Closed-loop simulation of an MPC applied to a cloth model
Original linear model by David Parent
Modified by Adria Luque so both COM and SOM are linear models,
and then added the nonlinear model to simulate reality
%}
clear; close all; clc;

addpath('..\required_files\cloth_model_New_NL\')
addpath('..\required_files\cloth_model_New_L\')
addpath('..\required_files\casadi-toolbox')
import casadi.*

plotAnim = 0;
animwWAM = 0;

% General Parameters
NTraj = 10;
nCOM = 4;
nSOM = 4;
nNLM = 10;
Hp = 25;
Ts = 0.020;
Wv = 0.2;
ExpSetN = 4;
NExp = 8;
NTrial = 2;
zsum0 = +0.01;
TCPOffset_local = [0; 0; 0.09];

% Opti parameters
xbound = 1;
ubound = 5*1e-3;
gbound = 0;  % 0 -> Equality constraint
W_Q = 1;
W_T = 1;
W_R = 10;

% Noise parameters
sigmaX = 0.075;

% -------------------


% Load trajectory to follow
phi_l_Traj = load(['../data/trajectories/ref_',num2str(NTraj),'L.csv']);
phi_r_Traj = load(['../data/trajectories/ref_',num2str(NTraj),'R.csv']);

% Get implied cloth size, position and angle wrt XZ
dphi_corners1 = phi_r_Traj(1,:) - phi_l_Traj(1,:);
lCloth = norm(dphi_corners1);
cCloth = (phi_r_Traj(1,:) + phi_l_Traj(1,:))/2 + [0 0 lCloth/2];
aCloth = atan2(dphi_corners1(2), dphi_corners1(1));

% Define COM parameters
COM = struct;
COM.row = nCOM;
COM.col = nCOM;
COM.mass = 0.1;
COM.grav = 9.8;
COM.dt = Ts;


% Load parameter table and select corresponding row
ThetaLUT = readtable('../learn_model/ThetaModelLUT.csv');

% Get the corresponding row(s)
LUT_COM_id = (ThetaLUT.ExpSetN == ExpSetN) & (ThetaLUT.NExp == NExp) & ...
             (ThetaLUT.NTrial == NTrial) & (ThetaLUT.Ts == Ts) & ...
             (ThetaLUT.nCOM == nCOM);
LUT_SOM_id = (ThetaLUT.NExp == NExp) & (ThetaLUT.NTrial == NTrial) & ... 
             (ThetaLUT.Ts == Ts) & (ThetaLUT.nCOM == nSOM);
LUT_COM = ThetaLUT(LUT_COM_id, :);
LUT_SOM = ThetaLUT(LUT_COM_id, :);

if (size(LUT_COM,1) > 1 || size(LUT_SOM,1) > 1)
    error("There are multiple rows with same experiment parameters.");
elseif (size(LUT_COM,1) < 1 || size(LUT_SOM,1) < 1)
    error("There are no saved experiments with those parameters.");
else
    thetaC = table2array(LUT_COM(:, contains(LUT_COM.Properties.VariableNames, 'Th_')));
    thetaS = table2array(LUT_SOM(:, contains(LUT_SOM.Properties.VariableNames, 'Th_')));
end

% Apply COM parameters
COM.stiffness = thetaC(1:3);
COM.damping = thetaC(4:6);
COM.z_sum = thetaC(7) + zsum0;


% Controlled coordinates (upper corners in x,y,z)
COM_node_ctrl = [nCOM*(nCOM-1)+1, nCOM^2];
COM.coord_ctrl = [COM_node_ctrl, ...
                  COM_node_ctrl+nCOM^2, ...
                  COM_node_ctrl+2*nCOM^2];

              
% Define the SOM (LINEAR)
SOM = struct;
SOM.row = nSOM;
SOM.col = nSOM;
SOM.mass = 0.1;
SOM.grav = 9.8;
SOM.dt = Ts;

% Real initial position in space
pos = create_lin_mesh(lCloth, nSOM, cCloth, aCloth);

% Apply SOM parameters
SOM.stiffness = thetaS(1:3);
SOM.damping = thetaS(4:6);
SOM.z_sum = thetaS(7) + zsum0;


% Controlled coordinates (upper corners in x,y,z)
SOM_node_ctrl = [nSOM*(nSOM-1)+1, nSOM^2];
SOM.coord_ctrl = [SOM_node_ctrl SOM_node_ctrl+nSOM^2 SOM_node_ctrl+2*nSOM^2];

% Define initial position of the nodes (needed for ext_force)
% Second half is velocity (initial v=0)
x_ini_SOM = [reshape(pos,[3*nSOM^2 1]); zeros(3*nSOM^2,1)];

% Reduce initial SOM position to COM size if necessary
[pos_rd,~] = take_reduced_mesh(x_ini_SOM(1:3*nSOM^2),x_ini_SOM(3*nSOM^2+1:6*nSOM^2), nSOM, nCOM);
x_ini_COM = [pos_rd; zeros(3*nCOM^2,1)];

% Rotate initial COM and SOM positions to XZ plane
RCloth_ini = [cos(aCloth) -sin(aCloth) 0; sin(aCloth) cos(aCloth) 0; 0 0 1];
posSOM_XZ = (RCloth_ini^-1 * pos')';
posCOM = reshape(x_ini_COM(1:3*nCOM^2), [nCOM^2,3]);
posCOM_XZ = (RCloth_ini^-1 * posCOM')';

% Initial position of the nodes
SOM.nodeInitial = lift_z(posSOM_XZ, SOM);
COM.nodeInitial = lift_z(posCOM_XZ, COM);


% Find initial spring length in each direction x,y,z
[SOM.mat_x, SOM.mat_y, SOM.mat_z] = compute_l0_linear(SOM,0);
[COM.mat_x, COM.mat_y, COM.mat_z] = compute_l0_linear(COM,0);

% Find linear matrices
[A_SOM, B_SOM, f_SOM] = create_model_linear_matrices(SOM);
[A_COM, B_COM, f_COM] = create_model_linear_matrices(COM);


% Third model as a real cloth representation (NL)
[NLM, pos_nl] = initialize_nl_model(lCloth,nNLM,cCloth,aCloth,Ts);
x_ini_NLM = [reshape(pos_nl,[3*nNLM^2 1]); zeros(3*nNLM^2,1)];

% Nonlinear model size = vision mesh size, must be reduced to SOM size
[pos_rd_nl,~] = take_reduced_mesh(x_ini_NLM(1:3*nNLM^2),x_ini_NLM(3*nNLM^2+1:6*nNLM^2), nNLM, nSOM);
x_ini_NLM_red = [pos_rd_nl; zeros(3*nSOM^2,1)];

% Nonlinear model corner coordinates
NLM_node_ctrl = [nNLM*(nNLM-1)+1, nNLM^2];
NLM.coord_ctrl = [NLM_node_ctrl NLM_node_ctrl+nNLM^2 NLM_node_ctrl+2*nNLM^2];
coord_lcNL = [1 nNLM 1+nNLM^2 nNLM^2+nNLM 2*nNLM^2+1 2*nNLM^2+nNLM]; 
n_states_nl = 3*2*nNLM^2;


%% Start casADi optimization problem
% Declare model variables
phi = SX.sym('phi',3*COM.row*COM.col);
dphi = SX.sym('dphi',3*COM.row*COM.col);
x = [phi; dphi];
n_states = length(x); % 3*2*nxC*nyC
u = SX.sym('u',6);

% Define model equations
x_next = SX.sym('xdot',6*COM.row*COM.col);
x_next(:) = A_COM*x + B_COM*u + COM.dt*f_COM;

% (x,u)->(x_next)
stfun = Function('stfun',{x,u},{x_next}); % nonlinear mapping function f(x,u)

% Lower corner coordinates for both models
coord_lcC = [1 nCOM 1+nCOM^2 nCOM^2+nCOM 2*nCOM^2+1 2*nCOM^2+nCOM]; 
coord_lcS = [1 nSOM 1+nSOM^2 nSOM^2+nSOM 2*nSOM^2+1 2*nSOM^2+nSOM];

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
ln = P([1 3 5],end) - Xk(coord_lcC([1,3,5])); %ln = left node
ln = abs(ln)./(norm(ln)+10^-6);
rn = P([2,4,6],end) - Xk(coord_lcC([2,4,6])); %rn = right node
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
    
    % Integrate until the end of the interval
    Xk_next = stfun(Xk,Uk);
    Xkn_r = Xk_next(coord_lcC);
    
    % Constant distance between upper corners
    Xkn_u = Xk_next(COM.coord_ctrl);
    g = [g; sum((Xkn_u([2,4,6]) - Xkn_u([1,3,5])).^2) - lCloth^2 ];
    lbg = [lbg; -gbound-1e-6];
    ubg = [ubg;  gbound+1e-6];
    
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
solver = nlpsol('solver', 'ipopt', nlp_prob,opts);


%----------------------------------%
%% MAIN SIMULATION LOOP EXECUTION %%
%----------------------------------%

% Initialize control
u_ini = x_ini_SOM(SOM.coord_ctrl);
u_bef = u_ini;
u_SOM = u_ini;

% Get Cloth orientation (rotation matrix)
cloth_x = u_SOM([2 4 6]) - u_SOM([1 3 5]);
cloth_y = [-cloth_x(2) cloth_x(1) 0]';
cloth_z = cross(cloth_x,cloth_y);

cloth_x = cloth_x/norm(cloth_x);
cloth_y = cloth_y/norm(cloth_y);
cloth_z = cloth_z/norm(cloth_z);
Rcloth = [cloth_x cloth_y cloth_z];
Rtcp = [cloth_y cloth_x -cloth_z];

% TCP initial position
tcp_ini = (u_SOM([1 3 5])+u_SOM([2 4 6]))'/2 + (Rcloth*TCPOffset_local)';

% Initialize things
reference = zeros(6*COM.row*COM.col, Hp+1);
store_somstate(:,1) = x_ini_SOM;
store_nlmstate(:,1) = x_ini_NLM;
store_u(:,1) = zeros(6,1);
store_pose(1) = struct('position', tcp_ini, ...
                       'orientation', rotm2quat(Rtcp));

tT0 = tic;
t0 = tic;
printX = 50;
for tk=2:size(phi_l_Traj,1)
    
    % Get new feedback value (eq. to "Spin once")
    x_noise_nl = [normrnd(0,sigmaX^2,[n_states_nl/2,1]); zeros(n_states_nl/2,1)];
    x_prev_noisy_nl = store_nlmstate(:,tk-1) + x_noise_nl;
    
    [phi_noisy, dphi_noisy] = take_reduced_mesh(x_prev_noisy_nl(1:3*nNLM^2), ...
                                      x_prev_noisy_nl(3*nNLM^2+1:6*nNLM^2), ...
                                      nNLM, nSOM);
    x_prev_noisy = [phi_noisy; dphi_noisy];
    x_prev_wavg = x_prev_noisy*Wv + store_somstate(:,tk-1)*(1-Wv);
    
    % The last Hp timesteps, trajectory should remain constant
    if tk>=size(phi_l_Traj,1)-Hp 
        Traj_l_Hp = repmat(phi_l_Traj(end,:), Hp,1);
        Traj_r_Hp = repmat(phi_r_Traj(end,:), Hp,1);
    else
        Traj_l_Hp = phi_l_Traj(tk:tk+Hp-1,:);
        Traj_r_Hp = phi_r_Traj(tk:tk+Hp-1,:);
    end
    
    % Rotate initial position to cloth base
    pos_ini_COM = reshape(x_ini_COM(1:3*nCOM^2),[nCOM^2,3]);
    vel_ini_COM = reshape(x_ini_COM(3*nCOM^2+1:6*nCOM^2),[nCOM^2,3]);
    
    pos_ini_COM_rot = (Rcloth^-1 * pos_ini_COM')';
    vel_ini_COM_rot = (Rcloth^-1 * vel_ini_COM')';
    x_ini_COM_rot = [reshape(pos_ini_COM_rot,[3*nCOM^2,1]);
                     reshape(vel_ini_COM_rot,[3*nCOM^2,1])];
                 
    % Rotate reference trajectory to cloth base
    Traj_l_Hp_rot = (Rcloth^-1 * Traj_l_Hp')';
    Traj_r_Hp_rot = (Rcloth^-1 * Traj_r_Hp')';
    
    % Define reference in the prediction horizon (sliding window)
    reference(:,1) = x_ini_COM_rot;
    reference(1,2:end) = Traj_l_Hp_rot(:,1)';
    reference(2,2:end) = Traj_r_Hp_rot(:,1)';
    reference(3,2:end) = Traj_l_Hp_rot(:,2)';
    reference(4,2:end) = Traj_r_Hp_rot(:,2)';
    reference(5,2:end) = Traj_l_Hp_rot(:,3)';
    reference(6,2:end) = Traj_r_Hp_rot(:,3)';
    
    % Initial seed of the optimization (for "u" and "r^a")
    % Add cloth size to z for LC->UC
    args_x0 = repmat([reference(1:6,end)+[0;0;0;0;lCloth;lCloth];zeros(6,1)],Hp,1);
    
    % Find the solution "sol"
    sol = solver('x0', args_x0, 'lbx', lbw, 'ubx', ubw, ...
                 'lbg', lbg, 'ubg', ubg, 'p', reference);

    % Get only controls from the solution
    % Control actions are upper corner displacements (incremental pos)
    % And they are in local base
    u_rot = reshape(full(sol.x(take_u)),6,Hp)'; 
    u_rot1 = u_rot(1,1:end)';
    u_rot2 = [u_rot1([1 3 5]) u_rot1([2 4 6])];
    
    % Convert back to global base
    u_lin2 = Rcloth * u_rot2;
    u_lin = reshape(u_lin2',[6,1]);
    
    % Output for Cartesian Ctrl is still u_SOM
    u_SOM = u_lin+u_bef;
    u_bef = u_SOM;
    
    % Linear SOM uses local variables too (rot)
    pos_ini_SOM = reshape(x_prev_wavg(1:3*nSOM^2), [nSOM^2,3]);
    vel_ini_SOM = reshape(x_prev_wavg(3*nSOM^2+1:6*nSOM^2), [nSOM^2,3]);
    pos_ini_SOM_rot = (Rcloth^-1 * pos_ini_SOM')';
    vel_ini_SOM_rot = (Rcloth^-1 * vel_ini_SOM')';
    x_ini_SOM_rot = [reshape(pos_ini_SOM_rot,[3*nSOM^2,1]);
                     reshape(vel_ini_SOM_rot,[3*nSOM^2,1])];
    
    % Simulate a SOM step
    next_state_SOM = A_SOM*x_ini_SOM_rot + B_SOM*u_rot1 + SOM.dt*f_SOM;
    
    % Simulate a NLM step
    [pos_nxt_NLM, vel_nxt_NLM] = simulate_cloth_step(x_prev_noisy_nl,u_SOM,NLM); 
    
    % Convert back to global axis
    pos_nxt_SOM_rot = reshape(next_state_SOM(1:3*nSOM^2), [nSOM^2,3]);
    vel_nxt_SOM_rot = reshape(next_state_SOM((1+3*nSOM^2):6*nSOM^2), [nSOM^2,3]); 
    pos_nxt_SOM = reshape((Rcloth * pos_nxt_SOM_rot')', [3*nSOM^2,1]);
    vel_nxt_SOM = reshape((Rcloth * vel_nxt_SOM_rot')', [3*nSOM^2,1]);
        
    % Close the loop
    [phired, dphired] = take_reduced_mesh(pos_nxt_SOM,vel_nxt_SOM, nSOM, nCOM);
    x_ini_COM = [phired; dphired];
    
    % Get new Cloth orientation (rotation matrix)
    cloth_x = u_SOM([2 4 6]) - u_SOM([1 3 5]);
    cloth_y = [-cloth_x(2) cloth_x(1) 0]';
    cloth_z = cross(cloth_x,cloth_y);
    
    cloth_x = cloth_x/norm(cloth_x);
    cloth_y = cloth_y/norm(cloth_y);
    cloth_z = cloth_z/norm(cloth_z);
    Rcloth = [cloth_x cloth_y cloth_z];
    Rtcp = [cloth_y cloth_x -cloth_z];
    
    % Real application with 1 robot: get EE pose
    TCPOffset = Rcloth * TCPOffset_local;
    PoseTCP = struct();
    PoseTCP.position = (u_SOM([1 3 5]) + u_SOM([2 4 6]))'/2 + TCPOffset';
    PoseTCP.orientation = rotm2quat(Rtcp);
    
    % Store things
    store_somstate(:,tk) = [pos_nxt_SOM; vel_nxt_SOM];
    store_nlmstate(:,tk) = [pos_nxt_NLM; vel_nxt_NLM];
    store_u(:,tk) = u_lin;
    store_pose(tk) = PoseTCP;
    
    if(mod(tk,printX)==0)
        t10 = toc(t0)*1000;
        fprintf(['Iter: ', num2str(tk), ...
            ' \t Avg. time/iter: ', num2str(t10/printX), ' ms \n']);
        t0 = tic;
    end
end
tT = toc(tT0);
time = 0:Ts:size(store_somstate,2)*Ts-Ts;
fprintf(['Total time: ',num2str(tT),' s \n']);


%% KPI
error_l = store_somstate(coord_lcS([1,3,5]),:)'-phi_l_Traj(1:end,1:3);
error_r = store_somstate(coord_lcS([2,4,6]),:)'-phi_r_Traj(1:end,1:3);

eMAE = mean(abs([error_l error_r]));
eMSE = mean([error_l error_r].^2);
eRMSE = sqrt(mean([error_l error_r].^2));

eMAEp  = mean([norm(eMAE([1,3,5]),2) norm(eMAE([2,4,6]),2)]);
eMSEp  = mean([norm(eMSE([1,3,5]),2) norm(eMSE([2,4,6]),2)]);
eRMSEp = mean([norm(eRMSE([1,3,5]),2) norm(eRMSE([2,4,6]),2)]);
eMAEm  = mean(eMAE,2); % Old "avg_error"
eMSEm  = mean(eMSE,2);
eRMSEm = mean(eRMSE,2);

% Save on struct
KPIs = struct();
KPIs.eMAE = eMAE;
KPIs.eMSE = eMSE;
KPIs.eRMSE = eRMSE;
KPIs.eMAEp = eMAEp;
KPIs.eMSEp = eMSEp;
KPIs.eRMSEp = eRMSEp;
KPIs.eMAEm = eMAEm;
KPIs.eMSEm = eMSEm;
KPIs.eRMSEm = eRMSEm;

% Display them
fprintf('\nExecution KPIs:\n');
fprintf(['- Coord. MAE:  \t', num2str(1000*eMAE),'\n']);
fprintf(['- Coord. RMSE: \t', num2str(1000*eRMSE),'\n']);
fprintf(['- Mean MAE:  \t', num2str(1000*eMAEm),' mm\n']);
fprintf(['- Norm MAE:  \t', num2str(1000*eMAEp),' mm\n']);
fprintf(['- Norm RMSE: \t', num2str(1000*eRMSEp),' mm\n']);



%% PLOT SOM CORNER EVOLUTION
fig1 = figure(1);
fig1.Color = [1,1,1];
fig1.Units = 'normalized';
fig1.Position = [0.5 0 0.5 0.9];

subplot(15,2,1:2:12);
plot(time, store_somstate(SOM.coord_ctrl([1 3 5]),:)','linewidth',1.5)
title('\textbf{Left upper corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(15,2,2:2:12);
plot(time, store_somstate(SOM.coord_ctrl([2 4 6]),:)','linewidth',1.5);
title('\textbf{Right upper corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(15,2,17:2:28);
plot(time, store_somstate(coord_lcS([1 3 5]),:)', 'linewidth',1.5);
hold on
plot(time, phi_l_Traj, '--k', 'linewidth',1.2);
hold off
title('\textbf{Left lower corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(15,2,18:2:28);
pa1som = plot(time, store_somstate(coord_lcS([2 4 6]),:)', 'linewidth',1.5);
hold on
pa1ref = plot(time, phi_r_Traj, '--k', 'linewidth',1.2);
hold off
title('\textbf{Right lower corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

Lgnd1 = legend([pa1som' pa1ref(1)], ...
               '$x_{SOM}$','$y_{SOM}$', '$z_{SOM}$', 'Ref', ...
               'Orientation','horizontal', 'Interpreter', 'latex');
Lgnd1.Position(1) = 0.5-Lgnd1.Position(3)/2;
Lgnd1.Position(2) = 0.06;


%% PLOT NLM CORNER EVOLUTION
fig2 = figure(2);
fig2.Color = [1,1,1];
fig2.Units = 'normalized';
fig2.Position = [0 0 0.5 0.9];

subplot(15,2,1:2:12);
plot(time, store_nlmstate(NLM.coord_ctrl([1 3 5]),:)','linewidth',1.5)
title('\textbf{Left upper corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(15,2,2:2:12);
plot(time, store_nlmstate(NLM.coord_ctrl([2 4 6]),:)','linewidth',1.5);
title('\textbf{Right upper corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(15,2,17:2:28);
plot(time, store_nlmstate(coord_lcNL([1 3 5]),:)', 'linewidth',1.5);
hold on
plot(time, phi_l_Traj, '--k', 'linewidth',1.2);
hold off
title('\textbf{Left lower corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(15,2,18:2:28);
pa2som = plot(time, store_nlmstate(coord_lcNL([2 4 6]),:)', 'linewidth',1.5);
hold on
pa1ref = plot(time, phi_r_Traj, '--k', 'linewidth',1.2);
hold off
title('\textbf{Right lower corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

Lgnd2 = legend([pa2som' pa1ref(1)], ...
               '$x_{NLM}$','$y_{NLM}$', '$z_{NLM}$', 'Ref', ...
               'Orientation','horizontal', 'Interpreter', 'latex');
Lgnd2.Position(1) = 0.5-Lgnd2.Position(3)/2;
Lgnd2.Position(2) = 0.06;


%% PLOT SOM IN 3D (W/ CLOTH MOVING)
fig3 = figure(3);
fig3.Color = [1,1,1];
fig3.Units = 'normalized';
fig3.Position = [0 0 0.5 0.90];

pov = [-30 20];
wampov = [-50 30];

SOM_ctrl = SOM.coord_ctrl(1:2);
SOM_lowc = coord_lcS(1:2);
store_sompos = store_somstate(1:3*nSOM^2,:);
All_uSOM = store_somstate(SOM.coord_ctrl,:);
TCP_pos = reshape([store_pose.position],[3,size(phi_l_Traj,1)])';
TCP_q   = reshape([store_pose.orientation],[4,size(phi_l_Traj,1)])';
TCP_rot = quat2rotm(TCP_q);
TCP_Tm  = [TCP_rot permute(TCP_pos',[1,3,2]);
          [0 0 0 1].*ones(1,1,size(phi_l_Traj,1))];
      
store_somx = store_sompos(1:nSOM^2,:);
limx = [floor(min(store_somx(:))*10), ceil(max(store_somx(:))*10)]/10;
store_somy = store_sompos(nSOM^2+1:2*nSOM^2,:);
limy = [floor(min(store_somy(:))*10), ceil(max(store_somy(:))*10)]/10;
store_somz = store_sompos(2*nSOM^2+1:3*nSOM^2,:);
limz = [floor(min(store_somz(:))*10), ceil(max(max(store_somz(:)), max(TCP_pos(:,3)))*10)]/10;

wamws = [-0.4 0.4 -0.8 0.2 -0.4 0.8];

plot3(All_uSOM(1:2,:)',All_uSOM(3:4,:)',All_uSOM(5:6,:)');
hold on
plot3(store_somx(SOM_lowc,:)', store_somy(SOM_lowc,:)', store_somz(SOM_lowc,:)');

scatter3(TCP_pos(1,1), TCP_pos(1,2), TCP_pos(1,3), 'om', 'filled');
plot3(TCP_pos(:,1), TCP_pos(:,2), TCP_pos(:,3), '--m');
plot3(phi_l_Traj(:,1),phi_l_Traj(:,2),phi_l_Traj(:,3), '--k');
plot3(phi_r_Traj(:,1),phi_r_Traj(:,2),phi_r_Traj(:,3), '--k');

scatter3(store_somx(:,1), store_somy(:,1), store_somz(:,1), '.b');
hold off
axis equal; box on; grid on;
xlim(limx);
ylim(limy);
zlim(limz);
set(gca, 'TickLabelInterpreter','latex');
xlabel('X', 'Interpreter','latex');
ylabel('Y', 'Interpreter','latex');
zlabel('Z', 'Interpreter','latex');
for fch=1:length(fig3.Children)
    if isa(fig3.Children(fch),'matlab.graphics.axis.Axes')
        fig3.Children(fch).View = pov;
    end
end

if(plotAnim > 0)
    if (animwWAM > 0)
        run("../with_robot/init_WAM.m");
    
        qini = wam.ikine(TCP_Tm(:,:,1), 'q0',qref);
        qt = qini;
    end
    
    pause(1);
    for t=2:size(store_somstate,2)

        scatter3(store_somx(:,t), store_somy(:,t), store_somz(:,t), '.b');
        hold on
        scatter3(store_somx(SOM_lowc,t), store_somy(SOM_lowc,t), ...
                 store_somz(SOM_lowc,t), 'ob');
        scatter3(TCP_pos(t,1), TCP_pos(t,2), TCP_pos(t,3), 'om', 'filled');
                
        plot3(TCP_pos(:,1), TCP_pos(:,2), TCP_pos(:,3), '--m');
        plot3(phi_l_Traj(:,1),phi_l_Traj(:,2),phi_l_Traj(:,3), '--k');
        plot3(phi_r_Traj(:,1),phi_r_Traj(:,2),phi_r_Traj(:,3), '--k');
        
        hold off
        axis equal; box on; grid on;
        xlim(limx);
        ylim(limy);
        zlim(limz);
        set(gca, 'TickLabelInterpreter','latex');
        xlabel('X', 'Interpreter','latex');
        ylabel('Y', 'Interpreter','latex');
        zlabel('Z', 'Interpreter','latex');
        
        if (animwWAM > 0)
            WAMbaseC = [0.8 0.8 0.8];
            hold on
            fill3([-0.1 0.1 0.1 -0.1],[-0.1 -0.1 0.1 0.1],[0 0 0 0],WAMbaseC,'FaceAlpha',0.7)
            fill3([-0.1 -0.1 -0.1 -0.1],[-0.1 -0.1 0.1 0.1],[-0.5 0 0 -0.5],WAMbaseC,'FaceAlpha',0.7)
            fill3([-0.1 -0.1 0.1 0.1],[-0.1 -0.1 -0.1 -0.1],[-0.5 0 0 -0.5],WAMbaseC,'FaceAlpha',0.7)
            fill3([0.1 0.1 0.1 0.1],[-0.1 -0.1 0.1 0.1],[-0.5 0 0 -0.5],WAMbaseC,'FaceAlpha',0.7)
            fill3([-0.1 -0.1 0.1 0.1],[0.1 0.1 0.1 0.1],[-0.5 0 0 -0.5],WAMbaseC,'FaceAlpha',0.7)
            hold off

            qt = wam.ikine(TCP_Tm(:,:,t), 'q0',qt);
            wam.plot(qt, 'workspace', wamws, ...
                         'notiles', 'noshadow', 'nobase', ...
                         'jointdiam', 0.6, 'jointlen', 0.8, ...
                         'lightpos', [-0.5 -0.5 1], 'fps', 30, ...
                         'linkcolor', [1 0.6 0], 'view', wampov, ...
                         'jointcolor', [0.7 0 1], 'pjointcolor', [0.7 0 1]);
        else
            for fch=1:length(fig3.Children)
                if isa(fig3.Children(fch),'matlab.graphics.axis.Axes')
                    fig3.Children(fch).View = pov;
                end
            end
            pause(1e-6);
        end
    end
    
    plot3(All_uSOM(1:2,:)',All_uSOM(3:4,:)',All_uSOM(5:6,:)');
    hold on
    plot3(store_somx(SOM_lowc,:)', store_somy(SOM_lowc,:)', store_somz(SOM_lowc,:)');

    scatter3(TCP_pos(t,1), TCP_pos(t,2), TCP_pos(t,3), 'om', 'filled');
    plot3(TCP_pos(:,1), TCP_pos(:,2), TCP_pos(:,3), '--m');
    plot3(phi_l_Traj(:,1),phi_l_Traj(:,2),phi_l_Traj(:,3), '--k');
    plot3(phi_r_Traj(:,1),phi_r_Traj(:,2),phi_r_Traj(:,3), '--k');

    scatter3(store_somx(:,t), store_somy(:,t), store_somz(:,t), '.b');
    hold off
    axis equal; box on; grid on;
    xlim(limx);
    ylim(limy);
    zlim(limz);
    set(gca, 'TickLabelInterpreter','latex');
    xlabel('X', 'Interpreter','latex');
    ylabel('Y', 'Interpreter','latex');
    zlabel('Z', 'Interpreter','latex');
    for fch=1:length(fig3.Children)
        if isa(fig3.Children(fch),'matlab.graphics.axis.Axes')
            fig3.Children(fch).View = pov;
        end
    end
    
end


%% PLOT NLM IN 3D
fig4 = figure(4);
fig4.Color = [1,1,1];
fig4.Units = 'normalized';
fig4.Position = [0.5 0 0.5 0.90];

NLM_ctrl = NLM.coord_ctrl(1:2);
NLM_lowc = coord_lcNL(1:2);
store_nlmpos = store_nlmstate(1:3*nNLM^2,:);
All_uNLM = store_nlmstate(NLM.coord_ctrl,:);

store_nlmx = store_nlmpos(1:nNLM^2,:);
store_nlmy = store_nlmpos(nNLM^2+1:2*nNLM^2,:);
store_nlmz = store_nlmpos(2*nNLM^2+1:3*nNLM^2,:);

plot3(All_uNLM(1:2,:)',All_uNLM(3:4,:)',All_uNLM(5:6,:)');
hold on
plot3(store_nlmx(NLM_lowc,:)', store_nlmy(NLM_lowc,:)', store_nlmz(NLM_lowc,:)');

scatter3(TCP_pos(1,1), TCP_pos(1,2), TCP_pos(1,3), 'om', 'filled');
plot3(TCP_pos(:,1), TCP_pos(:,2), TCP_pos(:,3), '--m');
plot3(phi_l_Traj(:,1),phi_l_Traj(:,2),phi_l_Traj(:,3), '--k');
plot3(phi_r_Traj(:,1),phi_r_Traj(:,2),phi_r_Traj(:,3), '--k');

scatter3(store_nlmx(:,1), store_nlmy(:,1), store_nlmz(:,1), '.b');
hold off
axis equal; box on; grid on;
xlim(limx);
ylim(limy);
zlim(limz);
set(gca, 'TickLabelInterpreter','latex');
xlabel('X', 'Interpreter','latex');
ylabel('Y', 'Interpreter','latex');
zlabel('Z', 'Interpreter','latex');
fig4.Children.View = pov;




