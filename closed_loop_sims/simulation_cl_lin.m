%{
Closed-loop simulation of an MPC applied to a linear cloth model
- Original linear model by David Parent, modified by Adrià Luque
- MPC and simulation by Adrià Luque
- Last Updated: September 2021
%}
clear; close all; clc;

% General Parameters
NTraj = 10;
Ts = 0.020;
Hp = 25;
nSOM = 4;
nCOM = 4;
TCPOffset_local = [0; 0; 0.09];

% Opti parameters
ubound = 50*1e-3;
gbound = 0; % (Eq. Constraint)
W_Q = 1;
W_R = 0.2;
opt_du  = 1;
opt_Qa  = 0;
opt_sto = 0;
opt_noise = 0;

% Noise parameters
sigmaD = opt_noise*0.003; % m/s
sigmaN = opt_noise*0.004; % m

% Plotting options
plotAnim = 0;
animwWAM = 0;
% ---------------------


% Add required directories, import CasADi
addpath('../required_files/cloth_model_New_L')
if (ispc)
    addpath('../required_files/casadi-toolbox-windows')
elseif (ismac)
    addpath('../required_files/casadi-toolbox-mac')
end
import casadi.*


% Load trajectory to follow
Ref_l = load(['../data/trajectories/ref_',num2str(NTraj),'L.csv']);
Ref_r = load(['../data/trajectories/ref_',num2str(NTraj),'R.csv']);
nPtRef = size(Ref_l,1);

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
    thetaC = table2array(LUT_COM(:, contains(LUT_COM.Properties.VariableNames,...
                                             'Th_')));
    thetaS = table2array(LUT_SOM(:, contains(LUT_SOM.Properties.VariableNames,...
                                             'Th_')));
end


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
COM.stiffness = thetaC(1:3);
COM.damping = thetaC(4:6);
COM.z_sum = thetaC(7);

% Important Coordinates (upper and lower corners in x,y,z)
COM_nd_ctrl = [nxC*(nyC-1)+1, nxC*nyC];
COM.coord_ctrl = [COM_nd_ctrl, COM_nd_ctrl+nxC*nyC, COM_nd_ctrl+2*nxC*nyC];
C_coord_lc = [1 nyC 1+nxC*nyC nxC*nyC+nyC 2*nxC*nyC+1 2*nxC*nyC+nyC]; 
COM.coord_lc = C_coord_lc;

              
% Define the SOM (LINEAR)
nxS = nSOM;
nyS = nSOM;
SOMlength = nxS*nyS;
SOM = struct;
SOM.row = nxS;
SOM.col = nyS;
SOM.mass = 0.1;
SOM.grav = 9.8;
SOM.dt = Ts;
SOM.stiffness = thetaS(1:3);
SOM.damping = thetaS(4:6);
SOM.z_sum = thetaS(7);

% Important Coordinates (upper and lower corners in x,y,z)
SOM_nd_ctrl = [nxS*(nyS-1)+1, nxS*nyS];
SOM.coord_ctrl = [SOM_nd_ctrl SOM_nd_ctrl+nxS*nyS SOM_nd_ctrl+2*nxS*nyS];
S_coord_lc = [1 nxS 1+nxS*nyS nxS*nyS+nxS 2*nxS*nyS+1 2*nxS*nyS+nxS];
SOM.coord_lc = S_coord_lc;

% Real initial position in space
pos = create_lin_mesh(lCloth, nSOM, cCloth, aCloth);

% Define initial position of the nodes (needed for ext_force)
% Second half is velocity (initial v=0)
x_ini_SOM = [reshape(pos,[3*nxS*nyS 1]); zeros(3*nxS*nyS,1)];

% Reduce initial SOM position to COM size if necessary
[reduced_pos,~] = take_reduced_mesh(x_ini_SOM(1:3*nxS*nyS), ...
                                    x_ini_SOM(3*nxS*nyS+1:6*nxS*nyS), ...
                                    nSOM, nCOM);
x_ini_COM = [reduced_pos; zeros(3*nxC*nyC,1)];

% Rotate initial COM and SOM positions to XZ plane
RCloth_ini = [cos(aCloth) -sin(aCloth) 0; sin(aCloth) cos(aCloth) 0; 0 0 1];
posSOM_XZ = (RCloth_ini^-1 * pos')';
posCOM = reshape(x_ini_COM(1:3*nxC*nyC), [nxC*nyC,3]);
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

Bd_COM = [zeros(3*nxC*nyC,3);
          [1 0 0].*ones(nxC*nyC,1);
          [0 1 0].*ones(nxC*nyC,1);
          [0 0 1].*ones(nxC*nyC,1)];
Bd_COM(3*nxC*nyC+COM.coord_ctrl,:) = 0;
Bd_SOM = [zeros(3*nxS*nyS,3);
          [1 0 0].*ones(nxS*nyS,1);
          [0 1 0].*ones(nxS*nyS,1);
          [0 0 1].*ones(nxS*nyS,1)];
Bd_SOM(3*nxS*nyS+SOM.coord_ctrl,:) = 0;


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
    lc_dist = Rp(:,end) - x0(COM.coord_lc);
    lc_dist = abs(lc_dist)/(norm(lc_dist)+eps);
    Q = diag(lc_dist);
end

for k = 1:Hp
    
    % Model Dynamics Constraint -> Definition
    x(:,k+1) = A_COM*x(:,k) + B_COM*u(:,k) + COM.dt*f_COM + ...
                              opt_sto*Bd_COM*d_hat(:,k);
    
    % Constraint: Constant distance between upper corners
    x_ctrl = x(COM.coord_ctrl,k+1);
    g = [g; sum((x_ctrl([2,4,6]) - x_ctrl([1,3,5])).^2) - lCloth^2 ];
    lbg = [lbg; -gbound];
    ubg = [ubg;  gbound];
    
    
    % Objective function
    x_err = x(COM.coord_lc,k+1) - Rp(:,k+1);
    objfun = objfun + x_err'*W_Q*Q*x_err;
    if (opt_du==0)
        objfun = objfun + u(:,k)'*W_R*u(:,k);
    else
        objfun = objfun + delta_u(:,k)'*W_R*delta_u(:,k);
    end
end

% Encapsulate in controller object
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
         'Ts = ',num2str(Ts*1000),' ms \t\t Hp = ',num2str(Hp),'\n', ...
         'nSOM = ',num2str(nSOM),' \t\t nCOM = ',num2str(nCOM),'\n', ...
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

% Initialize storage
in_params = zeros(2+6+3, max(n_states, Hp+1));
store_state(:,1) = x_ini_SOM;
store_noisy(:,1) = x_ini_SOM;
store_u(:,1) = zeros(6,1);
store_pose(1) = struct('position', tcp_ini, ...
                       'orientation', rotm2quat(Rtcp));

% Start timers and main loop
tT0 = tic;
t0 = tic;
printX = 100;
for tk=2:nPtRef
    
    % The last Hp+1 timesteps, trajectory should remain constant
    if tk>=nPtRef-(Hp+1) 
        Ref_l_Hp = repmat(Ref_l(end,:), Hp+1,1);
        Ref_r_Hp = repmat(Ref_r(end,:), Hp+1,1);
    else
        Ref_l_Hp = Ref_l(tk:tk+Hp,:);
        Ref_r_Hp = Ref_r(tk:tk+Hp,:);
    end
    
    % Rotate initial position to cloth base (R^T = R^-1)
    pos_ini_COM = reshape(x_ini_COM(1:3*nxC*nyC),[nxC*nyC,3]);
    vel_ini_COM = reshape(x_ini_COM(3*nxC*nyC+1:6*nxC*nyC),[nxC*nyC,3]);
    
    pos_ini_COM_rot = (Rcloth' * pos_ini_COM')';
    vel_ini_COM_rot = (Rcloth' * vel_ini_COM')';
    x_ini_COM_rot = [reshape(pos_ini_COM_rot,[3*nxC*nyC,1]);
                     reshape(vel_ini_COM_rot,[3*nxC*nyC,1])];
                 
    % Rotate reference trajectory to cloth base
    Ref_l_Hp_rot = (Rcloth' * Ref_l_Hp')';
    Ref_r_Hp_rot = (Rcloth' * Ref_r_Hp')';
    
    % Define input parameters for the optimizer (sliding window)
    in_params(1,:) = x_ini_COM_rot';
    in_params(2,1:6) = u_rot1';
    in_params(2+[1,3,5],1:Hp+1) = Ref_l_Hp_rot';
    in_params(2+[2,4,6],1:Hp+1) = Ref_r_Hp_rot';
    in_params(2+6+(1:3),1:Hp) = normrnd(0,sigmaD,[3,Hp]); % d_hat
    
    % Initial guess for optimizer (u: increments, guess UC=LC=Ref)
    args_x0 = [reshape(diff(in_params(2+(1:6),1:Hp),1,2),6*(Hp-1),1);zeros(6,1)];
    
    % Find the solution "sol"
    sol = controller('x0', args_x0, 'lbx', lbw, 'ubx', ubw, ...
                     'lbg', lbg, 'ubg', ubg, 'p', in_params);

    % Get only controls from the solution
    % Control actions are upper corner displacements (incremental pos)
    % And they are in local base
    u_rot = reshape(full(sol.x),6,Hp)'; 
    u_rot1 = u_rot(1,:)';
    u_rot2 = [u_rot1([1 3 5]) u_rot1([2 4 6])];
    
    % Convert back to global base
    u_lin2 = Rcloth * u_rot2;
    u_lin = reshape(u_lin2',[6,1]);
    
    % Output for Cartesian Ctrl is still u_SOM
    u_SOM = u_lin+u_bef;
    u_bef = u_SOM;
    
    % Add disturbance to SOM positions
    x_dist = Bd_SOM*normrnd(0,sigmaD,[3,1]);
    
    % Linear SOM uses local variables too (rot)  
    pos_ini_SOM = reshape(store_state(1:3*nxS*nyS,tk-1), [nxS*nyS,3]);
    vel_ini_SOM = reshape(store_state(3*nxS*nyS+1:6*nxS*nyS,tk-1),[nxS*nyS,3]);
    pos_ini_SOM_rot = (Rcloth' * pos_ini_SOM')';
    vel_ini_SOM_rot = (Rcloth' * vel_ini_SOM')';
    x_ini_SOM_rot = [reshape(pos_ini_SOM_rot,[3*nxS*nyS,1]);
                     reshape(vel_ini_SOM_rot,[3*nxS*nyS,1])];
    
    % Simulate a step
    next_state_SOM = A_SOM*x_ini_SOM_rot + B_SOM*u_rot1 + SOM.dt*f_SOM + x_dist;
    
    % Convert back to global axis
    pos_nxt_SOM_rot = reshape(next_state_SOM(1:3*nxS*nyS), [nxS*nyS,3]);
    vel_nxt_SOM_rot = reshape(next_state_SOM(1+3*nxS*nyS:6*nxS*nyS), [nxS*nyS,3]); 
    pos_nxt_SOM = reshape((Rcloth * pos_nxt_SOM_rot')', [3*nxS*nyS,1]);
    vel_nxt_SOM = reshape((Rcloth * vel_nxt_SOM_rot')', [3*nxS*nyS,1]);
    
    % Add sensor noise to positions
    pos_noise = normrnd(0,sigmaN,[3*nxS*nyS,1]);
    pos_noisy = pos_nxt_SOM + pos_noise*(tk>10);
        
    % Get COM states from SOM (Close the loop)
    [phired, dphired] = take_reduced_mesh(pos_noisy,vel_nxt_SOM, nSOM, nCOM);
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
    PoseTCP.position = (u_SOM([1 3 5]) + u_SOM([2 4 6]))' / 2 + TCPOffset';
    PoseTCP.orientation = rotm2quat(Rtcp);
    
    % Store things
    store_state(:,tk) = [pos_nxt_SOM; vel_nxt_SOM];
    store_noisy(:,tk) = [pos_noisy; vel_nxt_SOM];
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
fprintf(['-----------------------------------------\n', ...
         ' - Total time:  \t',num2str(tT),' s \n', ...
         ' - Avg. t/iter: \t',num2str(tT/nPtRef*1000),' ms \n']);


%% KPI
error_l = store_state(S_coord_lc([1,3,5]),:)'-Ref_l;
error_r = store_state(S_coord_lc([2,4,6]),:)'-Ref_r;

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
fprintf([' - Mean MAE:  \t\t', num2str(1000*eMAEm),' mm\n']);
fprintf([' - Norm RMSE: \t\t', num2str(1000*eRMSEp),' mm\n']);



%% PLOT CORNERS
time = 0:Ts:nPtRef*Ts-Ts;
plot_evo = store_state; %store_noisy | store_state 

fig1 = figure(1);
fig1.Color = [1,1,1];
fig1.Units = 'normalized';
fig1.Position = [0.5 0 0.5 0.9];

subplot(15,2,1:2:12);
plot(time, plot_evo(SOM.coord_ctrl([1 3 5]),:)','linewidth',1.5)
title('\textbf{Left upper corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(15,2,2:2:12);
plot(time, plot_evo(SOM.coord_ctrl([2 4 6]),:)','linewidth',1.5);
title('\textbf{Right upper corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(15,2,17:2:28);
plot(time, plot_evo(S_coord_lc([1 3 5]),:)', 'linewidth',1.5);
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
pa1som = plot(time, plot_evo(S_coord_lc([2 4 6]),:)', 'linewidth',1.5);
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


%% PLOT SOM IN 3D (W/ CLOTH MOVING)
fig3 = figure(3);
fig3.Color = [1,1,1];
fig3.Units = 'normalized';
fig3.Position = [0 0 0.5 0.90];

pov = [-30 20];
wampov = [-50 30];

SOM_ctrl = SOM.coord_ctrl(1:2);
SOM_lowc = S_coord_lc(1:2);
store_pos = plot_evo(1:3*SOMlength,:);
All_uSOM = plot_evo(SOM.coord_ctrl,:);
TCP_pos = reshape([store_pose.position],[3,size(Ref_l,1)])';
TCP_q   = reshape([store_pose.orientation],[4,size(Ref_l,1)])';
TCP_rot = quat2rotm(TCP_q);
TCP_Tm  = [TCP_rot permute(TCP_pos',[1,3,2]);
          [0 0 0 1].*ones(1,1,size(Ref_l,1))];

store_x = store_pos(1:SOMlength,:);
limx = [floor(min(store_x(:))*10), ceil(max(store_x(:))*10)]/10;
store_y = store_pos(SOMlength+1:2*SOMlength,:);
limy = [floor(min(store_y(:))*10), ceil(max(store_y(:))*10)]/10;
store_z = store_pos(2*SOMlength+1:3*SOMlength,:);
limz = [floor(min(store_z(:))*10), ceil(max(max(store_z(:)), max(TCP_pos(:,3)))*10)]/10;

wamws = [-0.4 0.4 -0.8 0.2 -0.4 0.8];

plot3(All_uSOM(1:2,:)',All_uSOM(3:4,:)',All_uSOM(5:6,:)');
hold on
plot3(store_x(SOM_lowc,:)', store_y(SOM_lowc,:)', store_z(SOM_lowc,:)');

scatter3(TCP_pos(1,1), TCP_pos(1,2), TCP_pos(1,3), 'om', 'filled');
plot3(TCP_pos(:,1), TCP_pos(:,2), TCP_pos(:,3), '--m');
plot3(Ref_l(:,1),Ref_l(:,2),Ref_l(:,3), '--k');
plot3(Ref_r(:,1),Ref_r(:,2),Ref_r(:,3), '--k');

scatter3(store_x(:,1), store_y(:,1), store_z(:,1), '.b');
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
    for t=2:nPtRef

        scatter3(store_x(:,t), store_y(:,t), store_z(:,t), '.b');
        hold on
        scatter3(store_x(SOM_lowc,t), store_y(SOM_lowc,t), ...
                 store_z(SOM_lowc,t), 'ob');
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
    plot3(store_x(SOM_lowc,:)', store_y(SOM_lowc,:)', store_z(SOM_lowc,:)');

    scatter3(TCP_pos(t,1), TCP_pos(t,2), TCP_pos(t,3), 'om', 'filled');
    plot3(TCP_pos(:,1), TCP_pos(:,2), TCP_pos(:,3), '--m');
    plot3(Ref_l(:,1),Ref_l(:,2),Ref_l(:,3), '--k');
    plot3(Ref_r(:,1),Ref_r(:,2),Ref_r(:,3), '--k');

    scatter3(store_x(:,t), store_y(:,t), store_z(:,t), '.b');
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








