function [Rwd, AllData] = simulation_cl_rtm_theta(theta, opts)

addpath('..\required_files\cloth_model_New_NL\')
addpath('..\required_files\cloth_model_New_L\')
addpath('..\required_files\casadi-toolbox')
import casadi.*

% Simulation Parameters
if nargin < 2
    NTraj = 6;
    Ts = 0.020;
    Hp = 25;
    Wv = 0.2;
    sigmaD = 0;
    sigmaN = 0;
    nCOM = 4;
    nSOM = 4;
    nNLM = 10;
    paramsSOM = [-300 -10 -225  -4 -2.5 -4 0.03];
    paramsCOM = [-300 -10 -225  -4 -2.5 -4 0.03];
    ubound = 50e-3;
    gbound = 0;
    opt_du = 1;
    opt_Qa = 0;
    RwdType = 1; % 1=RMSE, 2=Tov, 3=RMSE+Tov
else
    NTraj     = opts.NTraj;
    Ts        = opts.Ts;
    Hp        = opts.Hp;
    Wv        = opts.Wv;
    sigmaD    = opts.sigmaD;
    sigmaN    = opts.sigmaN;
    nCOM      = opts.nCOM;
    nSOM      = opts.nSOM;
    nNLM      = opts.nNLM;
    paramsSOM = opts.paramsSOM;
    paramsCOM = opts.paramsCOM;  
    ubound    = opts.ubound;
    gbound    = opts.gbound;
    opt_du    = opts.opt_du;
    opt_Qa    = opts.opt_Qa;
    RwdType   = opts.opt_Rwd;
end
% -------------------

% Extract input parameters
W_Q = theta.Q;
W_R = theta.R;


% Load trajectory to follow
Ref_l = load(['../data/trajectories/ref_',num2str(NTraj),'L.csv']);
Ref_r = load(['../data/trajectories/ref_',num2str(NTraj),'R.csv']);
nPtRef = size(Ref_l,1);

% Get implied cloth size, position and angle wrt XZ
dphi_corners1 = Ref_r(1,:) - Ref_l(1,:);
lCloth = norm(dphi_corners1);
cCloth = (Ref_r(1,:) + Ref_l(1,:))/2 + [0 0 lCloth/2];
aCloth = atan2(dphi_corners1(2), dphi_corners1(1));


% Define COM parameters
COM = struct;
COM.row = nCOM;
COM.col = nCOM;
COM.mass = 0.1;
COM.grav = 9.8;
COM.dt = Ts;
COM.stiffness = paramsCOM(1:3);
COM.damping = paramsCOM(4:6);
COM.z_sum = paramsCOM(7);

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
SOM.stiffness = paramsSOM(1:3);
SOM.damping = paramsSOM(4:6);
SOM.z_sum = paramsSOM(7);

% Important Coordinates (upper and lower corners in x,y,z)
SOM_node_ctrl = [nSOM*(nSOM-1)+1, nSOM^2];
SOM.coord_ctrl = [SOM_node_ctrl SOM_node_ctrl+nSOM^2 SOM_node_ctrl+2*nSOM^2];
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


% Third model as a real cloth representation (NL)
[NLM, pos_nl] = initialize_nl_model(lCloth,nNLM,cCloth,aCloth,Ts);
x_ini_NLM = [reshape(pos_nl,[3*nNLM^2 1]); zeros(3*nNLM^2,1)];
NL_coord_lc = NLM.coord_lc; 
n_states_nl = 3*2*nNLM^2;


%% Start casADi optimization problem

% Declare model variables
x = [SX.sym('pos',3*nCOM^2,Hp+1);
     SX.sym('vel',3*nCOM^2,Hp+1)];
u = SX.sym('u',6,Hp);
n_states = size(x,1); % 3*2*nxC*nyC

% Initial parameters of the optimization problem
P  = SX.sym('P', 2+6, max(n_states, Hp+1)); 
x0 = P(1, :)';
u0 = P(2, 1:6)';
Rp = P(2+(1:6), 1:Hp+1);

x(:,1) = x0;
delta_u = [u(:,1) - 0*u0, diff(u,1,2)];

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
    x(:,k+1) = (A_COM*x(:,k) + B_COM*u(:,k) + COM.dt*f_COM);
    
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

% Initialize things
in_params = zeros(2+6, max(n_states, Hp+1));
store_somstate(:,1) = x_ini_SOM;
store_nlmstate(:,1) = x_ini_NLM;
store_nlmnoisy(:,1) = x_ini_NLM;
store_u(:,1) = zeros(6,1);

tT0 = tic;
t0 = tic;
printX = floor(nPtRef/5);
for tk=2:nPtRef
    
    % Get new feedback value (eq. to "Spin once")
    x_noise_nl = [normrnd(0,sigmaN^2,[n_states_nl/2,1]); zeros(n_states_nl/2,1)];
    x_noisy_nl = store_nlmstate(:,tk-1) + x_noise_nl*(tk>10);
    
    [phi_noisy, dphi_noisy] = take_reduced_mesh(x_noisy_nl(1:3*nNLM^2), ...
                                      x_noisy_nl(3*nNLM^2+1:6*nNLM^2), ...
                                      nNLM, nSOM);
    x_noisy = [phi_noisy; dphi_noisy];
    somst_wavg = x_noisy*Wv + store_somstate(:,tk-1)*(1-Wv);
    
    % The last Hp timesteps, trajectory should remain constant
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
    
    % Add disturbance to NLM positions
    x_dist = [normrnd(0,sigmaD^2,[n_states_nl/2,1]); zeros(n_states_nl/2,1)];
    x_distd = store_nlmstate(:,tk-1) + x_dist*(tk>10);
    
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
    
    % Store things
    store_somstate(:,tk) = [pos_nxt_SOM; vel_nxt_SOM];
    store_nlmstate(:,tk) = [pos_nxt_NLM; vel_nxt_NLM];
    store_nlmnoisy(:,tk) = x_noisy_nl;
    store_u(:,tk) = u_lin;
    
    if(mod(tk,printX)==0)
        t10 = toc(t0)*1000;
        fprintf(['Iter: ', num2str(tk), ...
            ' \t Avg. time/iter: ', num2str(t10/printX), ' ms \n']);
        t0 = tic;
    end
end
tT = toc(tT0);
fprintf([' -- Total time: \t',num2str(tT),' s \n', ...
         ' -- Avg. t/iter: \t',num2str(tT/nPtRef*1000),' ms \n']);


%% KPI and Reward
error_som_l = 1000*(store_somstate(S_coord_lc([1,3,5]),:)'-Ref_l);
error_som_r = 1000*(store_somstate(S_coord_lc([2,4,6]),:)'-Ref_r);
error_nlm_l = 1000*(store_nlmstate(NL_coord_lc([1,3,5]),:)'-Ref_l);
error_nlm_r = 1000*(store_nlmstate(NL_coord_lc([2,4,6]),:)'-Ref_r);

error_l = error_som_l;
error_r = error_som_r;

%eMAE = mean(abs([error_l error_r]));
%eMAEp  = mean([norm(eMAE([1,3,5]),2) norm(eMAE([2,4,6]),2)]);
eRMSE = sqrt(mean([error_l error_r].^2));
eRMSEp = mean([norm(eRMSE([1,3,5]),2) norm(eRMSE([2,4,6]),2)]);
eTov = max(tT/(nPtRef*Ts) - 1, 0);

eRwds = -[eRMSEp eTov eRMSEp+eTov];
Rwd = eRwds(RwdType);

fprintf([' -- Sim. Reward: \t', num2str(Rwd), '\n']);
if isnan(Rwd)
    Rwd = -inf;
end


%% SAVE DATA
AllData = struct();
AllData.xSOM = store_somstate;
AllData.uSOM = store_somstate(SOM.coord_ctrl,:);
AllData.xNLM = store_nlmstate;
AllData.ulin = store_u;
AllData.eRMSE = eRMSEp;
AllData.eTov  = tT/(nPtRef*Ts);


end

