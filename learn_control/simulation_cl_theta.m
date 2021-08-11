function [Rwd, AllSt] = simulation_cl_theta(theta, opts)

addpath('..\required_files\cloth_model_FColtraro')
addpath('..\required_files\cloth_model_DParent')
addpath('..\required_files\casadi-toolbox')
import casadi.*

% Simulation Parameters
if nargin < 2
    % General Parameters
    nCOM = 4;
    nSOM = 10;
    NTraj = 3;
    Hp = 30;
    Ts = 0.010;

    % Opti parameters
    xbound = 1;
    ubound = 1*1e-3;
    gbound = 0;  % 0 -> Equality constraint
    W_Q = 1;
    W_T = 1;
    W_R = 10;
else
    nCOM = opts.nCOM;
    nSOM = opts.nSOM;
    NTraj = opts.NTraj;
    Hp = opts.Hp;
    Ts = opts.Ts;
    xbound = opts.xbound;
    ubound = opts.ubound;
    gbound = opts.gbound;
    W_Q = opts.W_Q;
    W_T = opts.W_T;
    W_R = opts.W_R;    
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
[reduced_pos,~] = take_reduced_mesh(x_ini_SOM(1:3*nxS*nyS),x_ini_SOM(3*nxS*nyS+1:6*nxS*nyS), nSOM, nCOM);
x_ini_COM = [reduced_pos; zeros(3*nxC*nyC,1)];

% Initial position of the nodes
COM.nodeInitial = lift_z(reshape(x_ini_COM(1:3*nxC*nyC), [nxC*nyC,3]), COM);


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
ln = abs(ln)./(norm(ln)+1e-6);
rn = P([2,4,6],end) - Xk(coord_l([2,4,6])); %rn = right node
rn = abs(rn)./(norm(rn)+1e-6);
Q = diag([ln(1) rn(1) ln(2) rn(2) ln(3) rn(3)]);

take_u=[];
i=1;

for k = 1:Hp
    % Artificial reference
    r_a = SX.sym(['r_a_' num2str(k)],6);
    w = [w; r_a];
    
    lbw = [lbw; -xbound*ones(6,1)]; 
    ubw = [ubw;  xbound*ones(6,1)];
    i=i+6;
    
    % New NLP variable for the control
    Uk = SX.sym(['U_' num2str(k)],6);
    w = [w; Uk];
    
    % The bounds of the control input are very relevant!! If too big, the COM
    %  allows it but the SOM not (the COM is a spring that from a time step
    %  to the next one can be mathematically infinetely stretched). The
    %  bounds should be approximately the maximum displacement expected.
    lbw = [lbw; -ubound*ones(6,1)]; 
    ubw = [ubw;  ubound*ones(6,1)];

    take_u=[take_u;(i:(i+6-1))'];
    i=i+6;
    
    % Integrate till the end of the interval
    Xk_next = stfun(Xk,Uk);
    Xkn_r = Xk_next(coord_l);
    
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
opt_settings = struct;
opt_settings.print_time = 0;
opt_settings.ipopt.print_level = 0; %0 to print the minimum, 3 to print the maximum
opt_settings.ipopt.warm_start_init_point = 'yes'; %warm start

solver = nlpsol('solver', 'ipopt', nlp_prob,opt_settings);

%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP
%
%%
% THE SIMULATION LOOP SHOULD START FROM HERE
%----------------------------------------------

% Define trajectory to follow
phi_l_Traj = load(['../data/trajectories/phi_',num2str(NTraj),'L.csv']);
phi_r_Traj = load(['../data/trajectories/phi_',num2str(NTraj),'R.csv']);

u_ini = x_ini_SOM(SOM.coord_ctrl);
u_bef = u_ini;

% Initialize things
reference = zeros(6*COM.row*COM.col, Hp+1);
store_state(:,1) = x_ini_SOM;
store_state(:,2) = x_ini_SOM;
store_u(:,1) = zeros(6,1);
store_u(:,2) = zeros(6,1);

t0 = tic;
printX = 100;
for tk=3:size(phi_l_Traj,1)
    
    % The last Hp timesteps, trajectory should remain constant
    if tk>=size(phi_l_Traj,1)-Hp 
        Traj_l_Hp = repmat(phi_l_Traj(end,:), Hp,1);
        Traj_r_Hp = repmat(phi_r_Traj(end,:), Hp,1);
    else
        Traj_l_Hp = phi_l_Traj(tk:tk+Hp-1,:);
        Traj_r_Hp = phi_r_Traj(tk:tk+Hp-1,:);
    end
    
    % Define reference in the prediction horizon (moving window)
    reference(:,1) = x_ini_COM;
    reference(1,2:end) = Traj_l_Hp(:,1)';
    reference(2,2:end) = Traj_r_Hp(:,1)';
    reference(3,2:end) = Traj_l_Hp(:,2)';
    reference(4,2:end) = Traj_r_Hp(:,2)';
    reference(5,2:end) = Traj_l_Hp(:,3)';
    reference(6,2:end) = Traj_r_Hp(:,3)';
    
    % Initial seed of the optimization (for "u" and "r^a")
    args_x0  = repmat([reference(1:6,end)+[0;0;0;0;0.3;0.3];zeros(6,1)],Hp,1);
    
    % Find the solution "sol"
    sol = solver('x0', args_x0, 'lbx', lbw, 'ubx', ubw, ...
                 'lbg', lbg, 'ubg', ubg, 'p', reference);

    % Get only controls from the solution
    u = reshape(full(sol.x(take_u)),6,Hp)'; 
    u_lin = u(1,1:end)';
    u_SOM = u_lin + u_bef;
    u_bef = u_SOM;
    
    next_state_SOM = cloth_simulator_secondorder([store_state(:,tk-1);store_state(:,tk-2)],u_SOM,SOM);
    
    phi_ini_SOM = next_state_SOM(1:3*nxS*nyS);
    dphi_ini_SOM = next_state_SOM((1+3*nxS*nyS):6*nxS*nyS);    
        
    % Close the loop
    [phired, dphired] = take_reduced_mesh(phi_ini_SOM,dphi_ini_SOM, nSOM, nCOM);
    x_ini_COM = [phired; dphired];
    
    % Store things
    store_state(:,tk) = [phi_ini_SOM; dphi_ini_SOM];
    store_u(:,tk) = u_lin;
    
    
    if(mod(tk,printX)==0)
        t10 = toc(t0)*1000;
        fprintf(['Iter: ', num2str(tk), ...
            ' \t Avg. time/iter: ', num2str(t10/printX), ' ms \n']);
        t0 = tic;
    end
end


%% COMPARE MODELS

All_StV = store_state;
All_StVrd = zeros(6*nxC*nyC, size(store_state,2));
All_uSOM = store_state(SOM.coord_ctrl,:);
All_ulin = store_u;
for i=1:size(store_state,2)
    pos_SOMi = store_state(1:3*SOMlength,i);
    vel_SOMi = store_state((1+3*SOMlength):6*SOMlength,i);
    [pos_rdi, vel_rdi] = take_reduced_mesh(pos_SOMi,vel_SOMi, nSOM, nCOM);
    All_StVrd(:,i) = [pos_rdi; vel_rdi];
end

All_StC = zeros(size(All_StVrd));
All_StC(:,1) = All_StVrd(:,1);
StCOM = All_StC(:,1);

for i=2:size(store_state,2)
    StCOM = A_COM*StCOM + B_COM*All_ulin(:,i) + Ts*f_COM;
    All_StC(:,i) = StCOM;
end

AllSt = struct();
AllSt.SOM = All_StV;
AllSt.SOMrd = All_StVrd;
AllSt.COM = All_StC;
AllSt.uSOM = All_uSOM;
AllSt.ulin = All_ulin;

avg_lin_error = mean(abs(All_StVrd-All_StC),2);
avg_lin_error_pos = avg_lin_error(1:3*COMlength);
Rwd = -norm(avg_lin_error_pos, 1);

disp(' ');
disp(['Reward: ', num2str(Rwd)]);

%% LOWER CORNER ERROR

avg_ref_error_1 = 1000*[mean(abs(store_state(coord_nl(1),:)'-phi_l_Traj(1:end,1)))
                    mean(abs(store_state(coord_nl(3),:)'-phi_l_Traj(1:end,2)))
                    mean(abs(store_state(coord_nl(5),:)'-phi_l_Traj(1:end,3)))];
avg_ref_error_4 = 1000*[mean(abs(store_state(coord_nl(2),:)'-phi_r_Traj(1:end,1)))
                    mean(abs(store_state(coord_nl(4),:)'-phi_r_Traj(1:end,2)))
                    mean(abs(store_state(coord_nl(6),:)'-phi_r_Traj(1:end,3)))];
avg_ref_error = [norm(avg_ref_error_1,2), norm(avg_ref_error_4,2)];

disp(['Tracking errors: ', num2str(avg_ref_error)]);

end

