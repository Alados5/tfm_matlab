function [Rwd, AllData] = simulation_cl_nl_theta(theta, opts)

addpath('..\required_files\cloth_model_New_NL\')
addpath('..\required_files\cloth_model_New_L\')
addpath('..\required_files\casadi-toolbox')
import casadi.*

% Simulation Parameters
if nargin < 2
    NTraj = 6;
    Ts = 0.020;
    Hp = 25;
    sigmaX = 0;
    nCOM = 4;
    nSOM = 4;
    paramsCOM = [-300 -10 -225  -4 -2.5 -4 0.03];
    xbound = 1.5;
else
    NTraj = opts.NTraj;
    Ts = opts.Ts;
    Hp = opts.Hp;
    sigmaX = opts.sigmaX;
    nCOM = opts.nCOM;
    nSOM = opts.nSOM;
    paramsCOM = opts.paramsCOM;  
    xbound = opts.xbound;
end
% -------------------

% Extract input parameters
ubound = theta(1);
gbound = theta(2);
W_Q = theta(3);
W_T = theta(4);
W_R = theta(5);


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
nxC = nCOM;
nyC = nCOM;
COM = struct;
COM.row = nxC;
COM.col = nyC;
COM.mass = 0.1;
COM.grav = 9.8;
COM.dt = Ts;

% Apply COM parameters
COM.stiffness = paramsCOM(1:3);
COM.damping = paramsCOM(4:6);
COM.z_sum = paramsCOM(7);


% Controlled coordinates (upper corners in x,y,z)
COM_node_ctrl = [nxC*(nyC-1)+1, nxC*nyC];
COM.coord_ctrl = [COM_node_ctrl, ...
                  COM_node_ctrl+nxC*nyC, ...
                  COM_node_ctrl+2*nxC*nyC];

              
% Define the SOM (NONLINEAR)
nxS = nSOM;
nyS = nSOM;
[SOM, pos] = initialize_nl_model(lCloth,nSOM,cCloth,aCloth,Ts);


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


%% Start casADi optimization problem
% Declare model variables
phi = SX.sym('phi',3*COM.row*COM.col);
dphi = SX.sym('dphi',3*COM.row*COM.col);
x = [phi; dphi];
n_states = length(x); %should be 96 in a 4x4 model
u = SX.sym('u',6);

% Define model equations
x_next = SX.sym('xdot',6*COM.row*COM.col);
x_next(:) = A_COM*x + B_COM*u + COM.dt*f_COM;

% (x,u)->(x_next)
stfun = Function('stfun',{x,u},{x_next}); % nonlinear mapping function f(x,u)

% Lower corner coordinates for both models
coord_l  = [1 nyC 1+nxC*nyC nxC*nyC+nyC 2*nxC*nyC+1 2*nxC*nyC+nyC]; 
coord_nl = [1 nxS 1+nxS*nyS nxS*nyS+nxS 2*nxS*nyS+1 2*nxS*nyS+nxS];

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
ln = P([1 3 5],end) - Xk(coord_l([1,3,5])); %ln = left node
ln = abs(ln)./(norm(ln)+10^-6);
rn = P([2,4,6],end) - Xk(coord_l([2,4,6])); %rn = right node
rn = abs(rn)./(norm(rn)+10^-6);
Q = diag([ln(1) rn(1) ln(2) rn(2) ln(3) rn(3)]);

% Q = [xl xr yl yr zl zr]
%Q = 0.5*diag([1 1 1 1 1 1]);


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
    Xkn_r = Xk_next(coord_l);
    
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
opts.ipopt.warm_start_init_point = 'yes'; %warm start

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);


%----------------------------------%
%% MAIN SIMULATION LOOP EXECUTION %%
%----------------------------------%

% Initialize control
u_ini = x_ini_SOM(SOM.coord_ctrl);
u_bef = u_ini;
u_SOM = u_ini;

% Get initial Cloth orientation (rotation matrix)
cloth_x = u_SOM([2 4 6]) - u_SOM([1 3 5]);
cloth_y = [-cloth_x(2) cloth_x(1) 0]';
cloth_z = cross(cloth_x,cloth_y);

cloth_x = cloth_x/norm(cloth_x);
cloth_y = cloth_y/norm(cloth_y);
cloth_z = cloth_z/norm(cloth_z);
Rcloth = [cloth_x cloth_y cloth_z];

% Initialize storage
reference = zeros(6*COM.row*COM.col, Hp+1);
store_state(:,1) = x_ini_SOM;
store_u(:,1) = zeros(6,1);

tT0 = tic;
t0 = tic;
printX = floor(nPtRef/5);
for tk=2:nPtRef
    
    % The last Hp timesteps, trajectory should remain constant
    if tk>=nPtRef-Hp 
        Traj_l_Hp = repmat(Ref_l(end,:), Hp,1);
        Traj_r_Hp = repmat(Ref_r(end,:), Hp,1);
    else
        Traj_l_Hp = Ref_l(tk:tk+Hp-1,:);
        Traj_r_Hp = Ref_r(tk:tk+Hp-1,:);
    end
    
    % Rotate initial position to cloth base
    pos_ini_COM = reshape(x_ini_COM(1:3*nxC*nyC),[nxC*nyC,3]);
    vel_ini_COM = reshape(x_ini_COM(3*nxC*nyC+1:6*nxC*nyC),[nxC*nyC,3]);
    
    pos_ini_COM_rot = (Rcloth^-1 * pos_ini_COM')';
    vel_ini_COM_rot = (Rcloth^-1 * vel_ini_COM')';
    x_ini_COM_rot = [reshape(pos_ini_COM_rot,[3*nxC*nyC,1]);
                     reshape(vel_ini_COM_rot,[3*nxC*nyC,1])];
                 
    % Rotate reference trajectory to cloth base
    Traj_l_Hp_rot = (Rcloth^-1 * Traj_l_Hp')';
    Traj_r_Hp_rot = (Rcloth^-1 * Traj_r_Hp')';
    
    % Define reference in the prediction horizon (moving window)
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
    
    % Add previous position for absolute position
    u_SOM = u_lin + u_bef;
    u_bef = u_SOM;
    
    % Get new Cloth orientation (rotation matrix)
    cloth_x = u_SOM([2 4 6]) - u_SOM([1 3 5]);
    cloth_y = [-cloth_x(2) cloth_x(1) 0]';
    cloth_z = cross(cloth_x,cloth_y);
    
    cloth_x = cloth_x/norm(cloth_x);
    cloth_y = cloth_y/norm(cloth_y);
    cloth_z = cloth_z/norm(cloth_z);
    Rcloth = [cloth_x cloth_y cloth_z];
    
    % Simulate a step of the SOM
    x_noise = [normrnd(0,sigmaX^2,[3*nSOM^2,1]); zeros(3*nSOM^2,1)];
    x_prev_noisy = store_state(:,tk-1) + x_noise;
    [pos_nxt_SOM, vel_nxt_SOM] = simulate_cloth_step(x_prev_noisy,u_SOM,SOM); 
    
    % Close the loop (update COM)
    [phired, dphired] = take_reduced_mesh(pos_nxt_SOM,vel_nxt_SOM, nSOM, nCOM);
    x_ini_COM = [phired; dphired];
    
    % Store things
    store_state(:,tk) = [pos_nxt_SOM; vel_nxt_SOM];
    store_u(:,tk) = u_lin;
    
    % Display progress
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
error_l = 1000*(store_state(coord_nl([1,3,5]),:)'-Ref_l);
error_r = 1000*(store_state(coord_nl([2,4,6]),:)'-Ref_r);

eMAE = mean(abs([error_l error_r]));
eRMSE = sqrt(mean([error_l error_r].^2));
eMAEp  = mean([norm(eMAE([1,3,5]),2) norm(eMAE([2,4,6]),2)]);
eRMSEp = mean([norm(eRMSE([1,3,5]),2) norm(eRMSE([2,4,6]),2)]);
eRwds = -[eMAEp eRMSEp];

Rwd = eRwds(2); %RMSE
fprintf([' -- Sim. Reward: \t', num2str(Rwd), '\n']);
if isnan(Rwd)
    Rwd = -inf;
end


%% SAVE DATA
AllData = struct();
AllData.xSOM = store_state;
AllData.uSOM = store_state(SOM.coord_ctrl,:);
AllData.ulin = store_u;


end

