%{
Closed-loop simulation of an MPC applied to cloth models
- Original linear model by David Parent, modified by Adrià Luque
- New NL model by Franco Coltraro
- MPC and simulation by Adrià Luque
- Both COM and SOM are linear, a third nonlinear model simulates reality
%}
clear; close all; clc;

% General Parameters
NTraj = 6;
Ts = 0.010;
Hp = 15;
Wv = 0.3;
nSOM = 4;
nCOM = 4;
nNLM = 10;
TCPOffset_local = [0; 0; 0.09];

% Opti parameters
ubound = 50*1e-3;
gbound = 0; % (Eq. Constraint)
W_Q = 0.05;
W_R = 1.00;
opt_du  = 1;
opt_Qa  = 0;
opt_sto = 0;

% Noise parameters
sigmaD = 0.003; % m/s
sigmaN = 0.003; % m

% Plotting options
plotAnim = 1;
animwWAM = 1;
plot_nlm = 0;
% ---------------------


% Add required directories, import CasADi
addpath('../required_files/cloth_model_New_L')
addpath('../required_files/cloth_model_New_NL')
if (ispc)
    addpath('../required_files/casadi-toolbox-windows')
elseif (ismac)
    disp('Download CasADi for Mac and add its path!');
end
import casadi.*


% Load trajectory to follow
Ref_l = load(['../data/trajectories/ref_',num2str(NTraj),'L.csv']);
Ref_r = load(['../data/trajectories/ref_',num2str(NTraj),'R.csv']);
nPtRef = size(Ref_l,1);
time = 0:Ts:nPtRef*Ts-Ts;

% Get implied cloth size, position and angle wrt XZ
dphi_corners1 = Ref_r(1,:) - Ref_l(1,:);
lCloth = norm(dphi_corners1);
cCloth = (Ref_r(1,:) + Ref_l(1,:))/2 + [0 0 lCloth/2];
aCloth = atan2(dphi_corners1(2), dphi_corners1(1));


% Load parameter table and select corresponding row(s)
ThetaLUT = readtable('../learn_model/ThetaMdl_LUT.csv');
LUT_SOM_id = (ThetaLUT.Ts == Ts) & (ThetaLUT.MdlSz == nSOM);
LUT_COM_id = (ThetaLUT.Ts == Ts) & (ThetaLUT.MdlSz == nCOM);
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


% Define COM parameters
COM = struct;
COM.row = nCOM;
COM.col = nCOM;
COM.mass = 0.1;
COM.grav = 9.8;
COM.dt = Ts;
COM.stiffness = thetaC(1:3);
COM.damping = thetaC(4:6);
COM.z_sum = thetaC(7);

% Important Coordinates (upper and lower corners in x,y,z)
COM_nd_ctrl = [nCOM*(nCOM-1)+1, nCOM^2];
COM.coord_ctrl = [COM_nd_ctrl, COM_nd_ctrl+nCOM^2, COM_nd_ctrl+2*nCOM^2];
C_coord_lc = [1 nCOM 1+nCOM^2 nCOM^2+nCOM 2*nCOM^2+1 2*nCOM^2+nCOM]; 
COM.coord_lc = C_coord_lc;


% Define the SOM (LINEAR)
SOM = struct;
SOM.row = nSOM;
SOM.col = nSOM;
SOM.mass = 0.1;
SOM.grav = 9.8;
SOM.dt = Ts;
SOM.stiffness = thetaS(1:3);
SOM.damping = thetaS(4:6);
SOM.z_sum = thetaS(7);

% Important Coordinates (upper and lower corners in x,y,z)
SOM_nd_ctrl = [nSOM*(nSOM-1)+1, nSOM^2];
SOM.coord_ctrl = [SOM_nd_ctrl SOM_nd_ctrl+nSOM^2 SOM_nd_ctrl+2*nSOM^2];
S_coord_lc = [1 nSOM 1+nSOM^2 nSOM^2+nSOM 2*nSOM^2+1 2*nSOM^2+nSOM];
SOM.coord_lc = S_coord_lc;

% Real initial position in space
pos = create_lin_mesh(lCloth, nSOM, cCloth, aCloth);

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

Bd_COM = [zeros(3*nCOM^2,3);
          [1 0 0].*ones(nCOM^2,1);
          [0 1 0].*ones(nCOM^2,1);
          [0 0 1].*ones(nCOM^2,1)];
Bd_COM(3*nCOM^2+COM.coord_ctrl,:) = 0;
Bd_SOM = [zeros(3*nSOM^2,3);
          [1 0 0].*ones(nSOM^2,1);
          [0 1 0].*ones(nSOM^2,1);
          [0 0 1].*ones(nSOM^2,1)];
Bd_SOM(3*nSOM^2+SOM.coord_ctrl,:) = 0;


% Third model as a real cloth representation (NL)
[NLM, pos_nl] = initialize_nl_model(lCloth,nNLM,cCloth,aCloth,Ts);
x_ini_NLM = [reshape(pos_nl,[3*nNLM^2 1]); zeros(3*nNLM^2,1)];
NL_coord_lc = NLM.coord_lc; 
n_states_nl = 3*2*nNLM^2;

Bd_NLM = [zeros(3*nNLM^2,3)
          [1 0 0].*ones(nNLM^2,1);
          [0 1 0].*ones(nNLM^2,1);
          [0 0 1].*ones(nNLM^2,1)];
Bd_NLM(3*nNLM^2+NLM.coord_ctrl,:) = 0;


%% Start casADi optimization problem

% Declare model variables
x = [SX.sym('pos',3*nCOM^2,Hp+1);
     SX.sym('vel',3*nCOM^2,Hp+1)];
u = SX.sym('u',6,Hp);
n_states = size(x,1); % 3*2*nxC*nyC

% Initial parameters of the optimization problem
P  = SX.sym('P', 2+6+3, max(n_states, Hp+1)); 
x0 = P(1, :)';
u0 = P(2, 1:6)';
Rp = P(2+(1:6), 1:Hp+1);
d_hat = P(2+6+(1:3), 1:Hp);

x(:,1) = x0;
delta_u = [u(:,1) - u0, diff(u,1,2)];

% Optimization variables
w = u(:);
lbw = -ubound*ones(6*Hp,1);
ubw = +ubound*ones(6*Hp,1);

% Initialize optimization variables
objfun = 0;     % Objective function
g = [];         % Constraints
lbg = [];       % Lower bounds of g
ubg = [];       % Upper bounds of g


% Adaptive weight calculation
if (opt_Qa == 0)
    % Disabled: Weigh coordinates constantly only with W_Q
    Q = 1;
else
    % Enabled: From the current LCpos to the desired at the horizon
    lc_dist = Rp(:,end) - x0(C_coord_lc);
    lc_dist = abs(lc_dist)/(norm(lc_dist)+eps);
    Q = diag(lc_dist);
end


for k = 1:Hp

    % Model Dynamics Constraint -> Definition
    x(:,k+1) = A_COM*x(:,k) + B_COM*u(:,k) + COM.dt*f_COM + opt_sto*Bd_COM*d_hat(:,k);
    
    % Constraint: Constant distance between upper corners
    x_ctrl = x(COM.coord_ctrl,k+1);
    g = [g; sum((x_ctrl([2,4,6]) - x_ctrl([1,3,5])).^2) - lCloth^2 ];
    lbg = [lbg; -gbound];
    ubg = [ubg;  gbound];
    
    
    % Objective function
    objfun = objfun + (x(C_coord_lc,k+1)-Rp(:,k+1))'*W_Q*Q*(x(C_coord_lc,k+1)-Rp(:,k+1));
    if (opt_du==0)
        objfun = objfun + u(:,k)'*W_R*u(:,k);
    else
        objfun = objfun + delta_u(:,k)'*W_R*delta_u(:,k);
    end
end


opt_prob = struct('f', objfun, 'x', w, 'g', g, 'p', P);
opt_config = struct;
opt_config.print_time = 0;
opt_config.ipopt.print_level = 0; %0 min print - 3 max
opt_config.ipopt.warm_start_init_point = 'yes'; %warm start

controller = nlpsol('ctrl_sol', 'ipopt', opt_prob, opt_config);


%----------------------------------%
%% MAIN SIMULATION LOOP EXECUTION %%
%----------------------------------%

% Initial info
fprintf(['Executing Reference Trajectory: ',num2str(NTraj), ...
         ' (',num2str(nPtRef),' pts) \n', ...
         'Ts = ',num2str(Ts*1000),' ms \t\t Hp = ',num2str(Hp), ...
         '\t \t \t Wv = ', num2str(Wv*100), '%% \n', ...
         'nSOM = ',num2str(nSOM),' \t\t nCOM = ',num2str(nCOM), ...
         '\t \t \t nNLM = ', num2str(nNLM), '\n', ...
         'lCloth = ',num2str(lCloth),' m \t aCloth = ',num2str(aCloth), ...
         ' rad \t cCloth = [', num2str(cCloth(1)), ', ' ...
         num2str(cCloth(2)),', ',num2str(cCloth(3)),'] m \n', ...
         '---------------------------------------------------------------\n']);

% Initialize control
u_ini  = x_ini_SOM(SOM.coord_ctrl);
u_lin  = zeros(6,1);
u_rot1 = u_lin;
u_bef  = u_ini;
u_SOM  = u_ini;

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

% Simulate some NLM steps to stabilize the NL model
warning('off','MATLAB:nearlySingularMatrix');
lastwarn('','');
[p_ini_NLM, ~] = simulate_cloth_step(x_ini_NLM,u_SOM,NLM);
[~, warnID] = lastwarn;
while strcmp(warnID, 'MATLAB:nearlySingularMatrix')
    lastwarn('','');
    x_ini_NLM = [p_ini_NLM; zeros(3*nNLM^2,1)];
    [p_ini_NLM, ~] = simulate_cloth_step(x_ini_NLM,u_SOM,NLM);
    [~, warnID] = lastwarn;
end
warning('on','MATLAB:nearlySingularMatrix');

% Initialize storage
in_params = zeros(2+6+3, max(n_states, Hp+1));
store_somstate(:,1) = x_ini_SOM;
store_nlmstate(:,1) = x_ini_NLM;
store_nlmnoisy(:,1) = x_ini_NLM;
store_u(:,1) = zeros(6,1);
store_pose(1) = struct('position', tcp_ini, ...
                       'orientation', rotm2quat(Rtcp));

tT = 0;
t1 = tic;
printX = 100;
for tk=2:nPtRef
    t0 = tic;
    
    % Get new noisy feedback value (eq. to "Spin once")
    x_noise_nl = [normrnd(0,sigmaN,[n_states_nl/2,1]); zeros(n_states_nl/2,1)];
    x_noisy_nl = store_nlmstate(:,tk-1) + x_noise_nl*(tk>10);
    
    [phi_noisy, dphi_noisy] = take_reduced_mesh(x_noisy_nl(1:3*nNLM^2), ...
                                       x_noisy_nl(3*nNLM^2+1:6*nNLM^2), ...
                                       nNLM, nSOM);
    x_noisy = [phi_noisy; dphi_noisy];
    somst_wavg = x_noisy*Wv + store_somstate(:,tk-1)*(1-Wv);
    
    % The last Hp+1 timesteps, trajectory should remain constant
    if tk>=nPtRef-(Hp+1) 
        Ref_l_Hp = repmat(Ref_l(end,:), Hp+1,1);
        Ref_r_Hp = repmat(Ref_r(end,:), Hp+1,1);
    else
        Ref_l_Hp = Ref_l(tk:tk+Hp,:);
        Ref_r_Hp = Ref_r(tk:tk+Hp,:);
    end
    
    % Get COM states from SOM (Close the loop)
    [phired, dphired] = take_reduced_mesh(somst_wavg(1:n_states/2), ...
                                          somst_wavg(n_states/2+1:n_states), ...
                                          nSOM, nCOM);
    x_ini_COM = [phired; dphired];
    
    % Rotate initial position to cloth base
    pos_ini_COM = reshape(x_ini_COM(1:3*nCOM^2),[nCOM^2,3]);
    vel_ini_COM = reshape(x_ini_COM(3*nCOM^2+1:6*nCOM^2),[nCOM^2,3]);
    
    pos_ini_COM_rot = (Rcloth^-1 * pos_ini_COM')';
    vel_ini_COM_rot = (Rcloth^-1 * vel_ini_COM')';
    x_ini_COM_rot = [reshape(pos_ini_COM_rot,[3*nCOM^2,1]);
                     reshape(vel_ini_COM_rot,[3*nCOM^2,1])];
                 
    % Rotate reference trajectory to cloth base
    Ref_l_Hp_rot = (Rcloth^-1 * Ref_l_Hp')';
    Ref_r_Hp_rot = (Rcloth^-1 * Ref_r_Hp')';
    
    % Define input parameters for the optimizer (sliding window)
    in_params(1,:) = x_ini_COM_rot';
    in_params(2,1:6) = u_rot1';
    in_params(2+[1,3,5],1:Hp+1) = Ref_l_Hp_rot';
    in_params(2+[2,4,6],1:Hp+1) = Ref_r_Hp_rot';
    in_params(2+6+(1:3),1:Hp) = normrnd(0,sigmaD,[3,Hp]); % d_hat
    
    % Initial guess for optimizer (u: increments, guess UC=LC=Ref)
    args_x0 = [reshape(diff(in_params(2+(1:6),1:Hp),1,2),6*(Hp-1),1); zeros(6,1)];
    
    % Find the solution "sol"
    sol = controller('x0', args_x0, 'lbx', lbw, 'ubx', ubw, ...
                     'lbg', lbg, 'ubg', ubg, 'p', in_params);

    % Get only controls from the solution
    % Control actions are upper corner displacements (incremental pos)
    % And they are in local base
    u_rot = reshape(full(sol.x),6,Hp)'; 
    u_rot1 = u_rot(1,1:end)';
    u_rot2 = [u_rot1([1 3 5]) u_rot1([2 4 6])];
    
    % Convert back to global base
    u_lin2 = Rcloth * u_rot2;
    u_lin = reshape(u_lin2',[6,1]);
    
    % Output for Cartesian Ctrl is still u_SOM
    u_SOM = u_lin+u_bef;
    u_bef = u_SOM;
    
    % Linear SOM uses local variables too (rot)
    pos_ini_SOM = reshape(somst_wavg(1:3*nSOM^2), [nSOM^2,3]);
    vel_ini_SOM = reshape(somst_wavg(3*nSOM^2+1:6*nSOM^2), [nSOM^2,3]);
    pos_ini_SOM_rot = (Rcloth^-1 * pos_ini_SOM')';
    vel_ini_SOM_rot = (Rcloth^-1 * vel_ini_SOM')';
    x_ini_SOM_rot = [reshape(pos_ini_SOM_rot,[3*nSOM^2,1]);
                     reshape(vel_ini_SOM_rot,[3*nSOM^2,1])];
    
    % Simulate a SOM step
    next_state_SOM = A_SOM*x_ini_SOM_rot + B_SOM*u_rot1 + SOM.dt*f_SOM;
    tT=tT+toc(t0);
    
    % Add disturbance to NLM positions
    x_dist = Bd_NLM*normrnd(0,sigmaD,[3,1]);
    x_distd = store_nlmstate(:,tk-1) + x_dist;
    
    % Simulate a NLM step
    [pos_nxt_NLM, vel_nxt_NLM] = simulate_cloth_step(x_distd,u_SOM,NLM); 
    
    % Convert back to global axis
    pos_nxt_SOM_rot = reshape(next_state_SOM(1:3*nSOM^2), [nSOM^2,3]);
    vel_nxt_SOM_rot = reshape(next_state_SOM((1+3*nSOM^2):6*nSOM^2), [nSOM^2,3]); 
    pos_nxt_SOM = reshape((Rcloth * pos_nxt_SOM_rot')', [3*nSOM^2,1]);
    vel_nxt_SOM = reshape((Rcloth * vel_nxt_SOM_rot')', [3*nSOM^2,1]);
    
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
    store_nlmnoisy(:,tk) = x_noisy_nl;
    store_u(:,tk) = u_lin;
    store_pose(tk) = PoseTCP;
    
    if(mod(tk,printX)==0)
        fprintf(['Iter: ', num2str(tk), ...
            ' \t Avg. time/iter: ', num2str(tT/tk*1000), ' ms \n']);
    end
end
tT = tT + toc(t0);
tT1 = toc(t1);
fprintf(['-----------------------------------------\n', ...
         ' [Times without NLM sim: extra ',num2str(tT1-tT),' s] \n',...
         ' - Total time:  \t',num2str(tT),' s \n', ...
         ' - Avg. t/iter: \t',num2str(tT/nPtRef*1000),' ms  \n']);


%% KPI
error_l = store_somstate(S_coord_lc([1,3,5]),:)'-Ref_l;
error_r = store_somstate(S_coord_lc([2,4,6]),:)'-Ref_r;

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
%fprintf([' - Coord. RMSE: \t', num2str(1000*eRMSE),'\n']);
fprintf([' - Mean MAE:  \t\t', num2str(1000*eMAEm),' mm\n']);
fprintf([' - Norm RMSE: \t\t', num2str(1000*eRMSEp),' mm\n']);



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
plot(time, store_somstate(S_coord_lc([1 3 5]),:)', 'linewidth',1.5);
hold on
plot(time, Ref_l, '--k', 'linewidth',1.2);
hold off
title('\textbf{Left lower corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(15,2,18:2:28);
pa1som = plot(time, store_somstate(S_coord_lc([2 4 6]),:)', 'linewidth',1.5);
hold on
pa1ref = plot(time, Ref_r, '--k', 'linewidth',1.2);
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

plot_nlevo = store_nlmnoisy;
if (plot_nlm == 1)
    fig2 = figure(2);
    fig2.Color = [1,1,1];
    fig2.Units = 'normalized';
    fig2.Position = [0 0 0.5 0.9];

    subplot(15,2,1:2:12);
    plot(time, plot_nlevo(NLM.coord_ctrl([1 3 5]),:)','linewidth',1.5)
    title('\textbf{Left upper corner}', 'Interpreter', 'latex')
    grid on
    xlabel('Time [s]', 'Interpreter', 'latex')
    ylabel('Position [m]', 'Interpreter', 'latex')
    xlim([0 time(end)])
    set(gca, 'TickLabelInterpreter', 'latex');

    subplot(15,2,2:2:12);
    plot(time, plot_nlevo(NLM.coord_ctrl([2 4 6]),:)','linewidth',1.5);
    title('\textbf{Right upper corner}', 'Interpreter', 'latex')
    grid on
    xlabel('Time [s]', 'Interpreter', 'latex')
    ylabel('Position [m]', 'Interpreter', 'latex')
    xlim([0 time(end)])
    set(gca, 'TickLabelInterpreter', 'latex');

    subplot(15,2,17:2:28);
    plot(time, plot_nlevo(NL_coord_lc([1 3 5]),:)', 'linewidth',1.5);
    hold on
    plot(time, Ref_l, '--k', 'linewidth',1.2);
    hold off
    title('\textbf{Left lower corner}', 'Interpreter', 'latex')
    grid on
    xlabel('Time [s]', 'Interpreter', 'latex')
    ylabel('Position [m]', 'Interpreter', 'latex')
    xlim([0 time(end)])
    set(gca, 'TickLabelInterpreter', 'latex');

    subplot(15,2,18:2:28);
    pa2som = plot(time, plot_nlevo(NL_coord_lc([2 4 6]),:)', 'linewidth',1.5);
    hold on
    pa1ref = plot(time, Ref_r, '--k', 'linewidth',1.2);
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
end


%% PLOT SOM IN 3D (W/ CLOTH MOVING)
fig3 = figure(3);
fig3.Color = [1,1,1];
fig3.Units = 'normalized';
fig3.Position = [0 0 0.5 0.90];

pov = [-30 20];
wampov = [-50 30];

SOM_ctrl = SOM.coord_ctrl(1:2);
SOM_lowc = S_coord_lc(1:2);
store_sompos = store_somstate(1:3*nSOM^2,:);
All_uSOM = store_somstate(SOM.coord_ctrl,:);
TCP_pos = reshape([store_pose.position],[3,nPtRef])';
TCP_q   = reshape([store_pose.orientation],[4,nPtRef])';
TCP_rot = quat2rotm(TCP_q);
TCP_Tm  = [TCP_rot permute(TCP_pos',[1,3,2]);
          [0 0 0 1].*ones(1,1,nPtRef)];
      
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
plot3(Ref_l(:,1),Ref_l(:,2),Ref_l(:,3), '--k');
plot3(Ref_r(:,1),Ref_r(:,2),Ref_r(:,3), '--k');

scatter3(store_somx(:,1), store_somy(:,1), store_somz(:,1), '.b');
hold off
axis equal; box on; grid on;
xlim(limx);
ylim(limy);
zlim(limz);
set(gca, 'TickLabelInterpreter','latex');
xlabel('$X$ [m]', 'Interpreter','latex');
ylabel('$Y$ [m]', 'Interpreter','latex');
zlabel('$Z$ [m]', 'Interpreter','latex');
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
        plot3(Ref_l(:,1),Ref_l(:,2),Ref_l(:,3), '--k');
        plot3(Ref_r(:,1),Ref_r(:,2),Ref_r(:,3), '--k');
        
        hold off
        axis equal; box on; grid on;
        xlim(limx);
        ylim(limy);
        zlim(limz);
        set(gca, 'TickLabelInterpreter','latex');
        xlabel('$X$ [m]', 'Interpreter','latex');
        ylabel('$Y$ [m]', 'Interpreter','latex');
        zlabel('$Z$ [m]', 'Interpreter','latex');
        
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
    plot3(Ref_l(:,1),Ref_l(:,2),Ref_l(:,3), '--k');
    plot3(Ref_r(:,1),Ref_r(:,2),Ref_r(:,3), '--k');

    scatter3(store_somx(:,t), store_somy(:,t), store_somz(:,t), '.b');
    hold off
    axis equal; box on; grid on;
    xlim(limx);
    ylim(limy);
    zlim(limz);
    set(gca, 'TickLabelInterpreter','latex');
    xlabel('$X$ [m]', 'Interpreter','latex');
    ylabel('$Y$ [m]', 'Interpreter','latex');
    zlabel('$Z$ [m]', 'Interpreter','latex');
    for fch=1:length(fig3.Children)
        if isa(fig3.Children(fch),'matlab.graphics.axis.Axes')
            fig3.Children(fch).View = pov;
        end
    end
    
end


%% PLOT NLM IN 3D

if (plot_nlm == 1)
    fig4 = figure(4);
    fig4.Color = [1,1,1];
    fig4.Units = 'normalized';
    fig4.Position = [0.5 0 0.5 0.90];

    NLM_ctrl = NLM.coord_ctrl(1:2);
    NLM_lowc = NL_coord_lc(1:2);
    store_nlmpos = plot_nlevo(1:3*nNLM^2,:);
    All_uNLM = plot_nlevo(NLM.coord_ctrl,:);

    store_nlmx = store_nlmpos(1:nNLM^2,:);
    store_nlmy = store_nlmpos(nNLM^2+1:2*nNLM^2,:);
    store_nlmz = store_nlmpos(2*nNLM^2+1:3*nNLM^2,:);

    plot3(All_uNLM(1:2,:)',All_uNLM(3:4,:)',All_uNLM(5:6,:)');
    hold on
    plot3(store_nlmx(NLM_lowc,:)', store_nlmy(NLM_lowc,:)', store_nlmz(NLM_lowc,:)');

    scatter3(TCP_pos(1,1), TCP_pos(1,2), TCP_pos(1,3), 'om', 'filled');
    plot3(TCP_pos(:,1), TCP_pos(:,2), TCP_pos(:,3), '--m');
    plot3(Ref_l(:,1),Ref_l(:,2),Ref_l(:,3), '--k');
    plot3(Ref_r(:,1),Ref_r(:,2),Ref_r(:,3), '--k');

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
end




