%% STATEMENT
%{
1. Compute the DH parameters. Set numerical values to the link lengths and
offsets

2. With the aid of the Symbolic Toolbox:
  a) Implement a script/function to compute the Forward Kinematics
  b) Implement a script/function to compute the Geomatric and Analytic
  Jacobians
  c) Implement a script/function to compute the Inverse Kinematics using a
  closed form
  d) Implement a script/function to compute the Inverse Kinematics using an
  iterative form

3. Using the Robotics Toolbox, model the manipulator as a SerialLink
instance, and create a script file to verify the correctness of the
functions implemented in the previous question.


%}
clear; clc; close all;


%% 1. DH PARAMETERS

d1 = 0.75;
d2 = 0;
th3 = 0;
d4 = -0.15;
a1 = 0.4;
a2 = 0.3;
a3 = 0;
a4 = 0;
alpha1 = 0;
alpha2 = pi;
alpha3 = -pi;
alpha4 = pi;

DH = [  0  d1  a1  alpha1;
        0  d2  a2  alpha2;
      th3   0  a3  alpha3;
        0  d4  a4  alpha4];
  


%% 2. Symbolic Toolbox Model

%% 2a) Forward Kinematics

% Create model
[T04, q, Ai4] = getKinModel(DH, 'rrpr');

% Same commands on script for specific case:
%{
% Joint variables
syms qi q1 q2 q3 q4
assume([qi q1 q2 q3 q4],'real')

% General DH Parameters
syms ai di alphai
assume([ai di alphai],'real')

% General A matrix
Ai = [cos(qi) -sin(qi)*cos(alphai)  sin(qi)*sin(alphai) ai*cos(qi);
      sin(qi)  cos(qi)*cos(alphai) -cos(qi)*sin(alphai) ai*sin(qi);
         0     sin(alphai)          cos(alphai)         di;
         0     0                    0                   1];

% Transform matrices
A01 = subs(Ai, [qi,di,ai,alphai], [q1,  d1, a1, alpha1]);
A12 = subs(Ai, [qi,di,ai,alphai], [q2,  d2, a2, alpha2]);
A23 = subs(Ai, [qi,di,ai,alphai], [th3, q3, a3, alpha3]);
A34 = subs(Ai, [qi,di,ai,alphai], [q4,  d4, a4, alpha4]);

% Final FK
T04  = simplify(A01*A12*A23*A34);
q = [q1 q2 q3 q4];
%}

% Test one configuration
qI = [pi/2 -pi/2 0.4 pi/2];
[TI, xI] = scaraFK(T04, qI);



%% 2b) Jacobians

% 2b I) Geometric Jacobian

J = getJacob(Ai4);

% Same commands on script for specific case:
%{
% Origins and z vectors of the frames
p0 = [0 0 0 1]';
p1 = A01*p0;
p2 = A01*A12*p0;
p3 = A01*A12*A23*p0;
p4 = A01*A12*A23*A34*p0;
z0 = [0 0 1]';
z1 = A01(1:3,1:3)*z0;
z2 = A01(1:3,1:3)*A12(1:3,1:3)*z0;
z3 = A01(1:3,1:3)*A12(1:3,1:3)*A23(1:3,1:3)*z0;

% Linear part
J11 = cross(z0, p4(1:3)-p0(1:3));
J12 = cross(z1, p4(1:3)-p1(1:3));
J14 = cross(z3, p4(1:3)-p3(1:3));

% Geometric Jacobian
J = [J11, J12, z2 J14; 
      z0,  z1,  0  z3];
        
J = simplify(J);
%}

% Test for same configuration
JI = double(subs(J, q, qI));
JI2 = getJacob(double(subs(Ai4,q,qI)),'rrpr');


% 2b II) Analytic Jacobian
theta = pi;

phi1 = q(1) + q(2) + q(4);
psi1 = 0;
JAphi = getJacobA(T04, [phi1; theta; psi1]);

phi2 = 0;
psi2 = -q(1) - q(2) - q(4);
JApsi = getJacobA(T04, [phi2; theta; psi2]);

% Same commands on script for specific case:
%{
q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
% Linear/Angular separation
Pe = T04(1:3,4);
TPhie = T04(1:3, 1:3);

% Euler angles (ZYZ)
phi = q(1) + q(2) + q(4);
theta = pi;
psi = 0;
Phie = [phi;theta;psi];

% Linear part
dPedq1 = diff(Pe, q1);
dPedq2 = diff(Pe, q2);
dPedq3 = diff(Pe, q3);
dPedq4 = diff(Pe, q4);
Jp = [dPedq1 dPedq2 dPedq3 dPedq4];
Jp = simplify(Jp);

% Angular part
dPhiedq1 = diff(Phie, q1);
dPhiedq2 = diff(Phie, q2);
dPhiedq3 = diff(Phie, q3);
dPhiedq4 = diff(Phie, q4);
JPhi = [dPhiedq1 dPhiedq2 dPhiedq3 dPhiedq4];
JPhi = simplify(JPhi);

JA = [Jp; JPhi]; 

TPhi = [0 -sin(phi) cos(phi)*sin(theta);
        0  cos(phi) sin(phi)*sin(theta);
        1     0              cos(theta)];
%}

% ZYZ Euler angles -> TA (J = TA*JA check)
TPhi1 = [0 -sin(phi1)  0;
         0  cos(phi1)  0;
         1     0      -1];
     
TPhi2 = [0 -sin(phi2)  0;
         0  cos(phi2)  0;
         1     0      -1];

TA1 = [eye(3),  zeros(3);
     zeros(3),   TPhi1];
 
TA2 = [eye(3),  zeros(3);
     zeros(3),   TPhi2];
 
Jcheck1 = TA1*JAphi;
Jcheck2 = TA1*JApsi;

Jeq1 = isequal(J, Jcheck1);
Jeq2 = isequal(J, Jcheck2);


%% 2c) Inverse Kinematics (closed)

% Test same configuration
qI_IK = scaraIK(TI, DH);
TI_IK = scaraFK(T04, qI_IK);

eIK = max(max(abs(TI-TI_IK)));

% Test new functions with Euler angles
[~, xIphi] = scaraFK(T04, qI, 'phi');
[~, xIpsi] = scaraFK(T04, qI, 'psi');
qIgeo_IK = scaraIKx(xI,DH,-1);
qIphi_IK = scaraIKx(xIphi,DH,-1,'phi');
qIpsi_IK = scaraIKx(xIpsi,DH,-1,'psi');

eIK_geo = max(abs(qIgeo_IK-qI));
eIK_phi = max(abs(qIphi_IK-qI));
eIK_psi = max(abs(qIpsi_IK-qI));



%% 2d) Inverse Kinematics (iterative)

% Methods:
%   JA pseudoinverse: dq = pinv(JA(q))*k*e
%   JA transpose: dq = JA(q)'*k*e

% Desired position obtention
qd = [pi/2 -pi/3 0.3 -pi];
[Td, xd] = scaraFK(T04, qd, 'phi');

% Main parameters
q0 = [0 0 0 0];
K = 100;
eStop = 5e-4;
JA = JAphi;
Ts = 1e-3;
TfMax = 10;
Eul = 'phi';


% Implementations:
%   I)  Function
%   II) Simulink model


% 2d I) scaraIKiter function
DoIKiterfcn = logical(input('Call iterative IK function [1/0]? '));
if DoIKiterfcn
    JAMethod = logical(input('Fcn IK Method [0: pinv(JA) | 1: JA^T]? '));
    % Input parameters:     (xd, Tq, JAq, Ts, q0, JAMethod, K, eStop, Eul, Timeout)
    try qd_IKiterf = scaraIKiter(xd, T04, JA, Ts, q0, JAMethod, K, eStop, Eul, TfMax);
    clc;
    Td_IKiterf = scaraFK(T04, qd_IKiterf);
    eIKiterf = max(max(abs(Td-Td_IKiterf)));
    catch
        disp('Function did not converge.');
    end
end
disp(' ');


% 2d II) Simulink IK_iterative model
DoIKi = logical(input('Continue with IK simulation and RTB checks [1/0]? '));
if ~DoIKi, return; end

% Select method
NotSelected = 1;
while NotSelected
    try JAMethod = logical(input('Sim IK Method [0: pinv(JA) | 1: JA^T]? '));
        disp(['Selected: ',num2str(JAMethod)]);
        disp(' ');
        NotSelected=0;
    catch
        disp('Invalid input');
        disp(' ');
    end
end

% Eul conversion
if strcmp(Eul,'phi'), EulSim = 1; Thz = 4; Thtxt = '\varphi';
elseif strcmp(Eul, 'psi'), EulSim = 3; Thz = 6; Thtxt = '\psi';
else, EulSim = 0; Thz = 6; Thtxt = '\phi';
end

% Simulate algorithm
tic;
sim('IK_iterative.slx');
tsimf = toc;
Tf = qout.time(end);
disp(tsimf);
disp(' ');

xf_iter = xe.data(end,:)';
qf_iter = qout.data(end,:);
eIKiters = max(abs(xd-xf_iter));


% Plot results
DoPlots = logical(input('Plot iterative IK results [1/0]? '));
if DoPlots
    close all;
    hf1 = figure(1);
    hf1.Color = [1 1 1];
    hf1.Position = [000 700 600 600];
    set(hf1, 'menubar', 'none')
    hf1s1 = subplot(2,1,1);
    plot([0 Tf], xd(1)*ones(1:2), '--')
    hold on
    plot([0 Tf], xd(2)*ones(1:2), '--')
    plot([0 Tf], xd(3)*ones(1:2), '--')
    plot(0, 0, '--') % Placeholder for legend
    plot(xe.time, xe.data(:,1), 'b', 'LineWidth', 1);
    plot(xe.time, xe.data(:,2), 'r', 'LineWidth', 1);
    plot(xe.time, xe.data(:,3), 'Color', [1 0.6 0], 'LineWidth', 1);
    plot(0, 0, 'Color', [0.6 0 1]) % Placeholder for legend
    hold off
    grid on
    xlabel('Simulated time [s]','Interpreter','latex','FontSize',12);
    ylabel('EE Position [m]', 'Interpreter','latex','FontSize',12);
    title('\textbf{End Efector position and orientation: real ($x_e(t)$) and desired ($x_d(t)$)}', ...
          'Interpreter', 'latex','FontSize',14);
    legend('$x_{d1} = p_{xd} \quad$', '$x_{d2} = p_{yd} \quad$', '$x_{d3} = p_{zd} \quad$', ...
          ['$x_{d',num2str(Thz),'} = ',Thtxt,'_d \quad$'], ...
           '$x_{e1} = p_{xe} \quad$', '$x_{e2} = p_{ye} \quad$', '$x_{e3} = p_{ze} \quad$', ...
          ['$x_{e',num2str(Thz),'} = ',Thtxt,'_e \quad$'], ...
           'Interpreter', 'latex','FontSize',12, 'NumColumns', 4, ...
           'Location','SouthOutside', 'Orientation', 'horizontal');
       
    hf1s2 = subplot(2,1,2);
    plot([0 Tf], xd(Thz)*ones(1:2), '--', 'Color', [0.4940 0.1840 0.5560])
    hold on
    plot(xe.time, xe.data(:,Thz), 'Color', [0.6 0 1], 'LineWidth', 1);
    hold off
    grid on
    xlabel('Simulated time [s]','Interpreter','latex','FontSize',12);
    ylabel('EE Orientation [rad]', 'Interpreter','latex','FontSize',12);
    
    hf1s1.Position = [0.1 0.55 0.8 0.4];
    hf1s2.Position = [0.1 0.1 0.8 0.3];
    
    
    hf2 = figure(2);
    hf2.Color = [1 1 1];
    hf2.Position = [800 700 600 600];
    set(hf2, 'menubar', 'none')
    subplot(3,1,1:2)
    plot(xe.time, xe.data(:,1)-xd(1), 'b', 'LineWidth', 1);
    hold on
    plot(xe.time, xe.data(:,2)-xd(2), 'r', 'LineWidth', 1);
    plot(xe.time, xe.data(:,3)-xd(3), 'Color', [1 0.6 0], 'LineWidth', 1);
    plot(xe.time, xe.data(:,Thz)-xd(Thz), 'Color', [0.6 0 1], 'LineWidth', 1);
    hold off
    grid on
    xlabel('Simulated time [s]','Interpreter','latex','FontSize',12);
    ylabel('Error: [m] or [rad]', 'Interpreter','latex','FontSize',12);
    title('\textbf{End Efector position and orientation errors ($e(t)$)}', ...
          'Interpreter', 'latex','FontSize',14);
    legend('$e_{1} = e_{x}$', '$e_{2} = e_{y}$', ...
           '$e_{3} = e_{z}$', ['$e_{',num2str(Thz),'} = e_',Thtxt,'$'], ...
           'Interpreter', 'latex','FontSize',12);
    
    subplot(3,1,3)
    plot(xe.time, xe.data(:,1)-xd(1), 'b', 'LineWidth', 1);
    hold on
    plot(xe.time, xe.data(:,2)-xd(2), 'r', 'LineWidth', 1);
    plot(xe.time, xe.data(:,3)-xd(3), 'Color', [1 0.6 0], 'LineWidth', 1);
    plot(xe.time, xe.data(:,Thz)-xd(Thz), 'Color', [0.6 0 1], 'LineWidth', 1);
    hold off
    grid on
    ylim([-10*eStop 10*eStop]);
    xlabel('Simulated time [s]','Interpreter','latex','FontSize',12);
    ylabel('Error: [m] or [rad]', 'Interpreter','latex','FontSize',12);
    title('\textbf{Detail of the errors ($e(t)$)}', ...
          'Interpreter', 'latex','FontSize',14);
    
end

SavePlots = logical(input('Save plots as separate files [1/0]? '));
if SavePlots
    if isempty(dir('Plots')), mkdir('Plots'); end
    
    saveas(hf1,['Plots/FW_1_xEE_JAM', num2str(JAMethod),'.png']);
    saveas(hf2,['Plots/FW_2_xerror_JAM', num2str(JAMethod),'.png']);
    
end
disp(' ');



%% 3. Robotics Toolbox check

% Create SCARA manipulator
L(1) = Link([DH(1,:), 0]);
L(2) = Link([DH(2,:), 0]);
L(3) = Link([DH(3,:), 1]);
L(4) = Link([DH(4,:), 0]);
scara = SerialLink(L, 'name', 'SCARA [ALA]');
scara.links(3).qlim = [0 d1-d4];

% Test and compare a configuration
qr = qd; %qd = [pi/2 -pi/3 0.3 -pi];

[Tr, xr] = scaraFK(T04,qr);
Tr_rtb = scara.fkine(qr);
eTr = max(max(abs(Tr-double(Tr_rtb))));

Jr = double(subs(J, q, qr));
Jr_rtb  = scara.jacob0(qr);
eJr = max(max(abs(Jr-Jr_rtb)));

TPhi1r = double(subs(TPhi1,q,qr));
TPhi2r = double(subs(TPhi2,q,qr));
TPhi_rtb = eul2jac(Tr_rtb.toeul);
eTPhi = max(max(min(abs(TPhi_rtb - TPhi1r), abs(TPhi_rtb - TPhi2r))));

qrIK1 = scaraIK(Tr, DH, +1);
qrIK2 = scaraIK(Tr, DH, -1);
qrIK_rtb = scara.ikine(Tr_rtb, 'mask', [1 1 1 0 0 1], ...
                                 'q0', [1 -1 0 -3]);
eIKr = max(min(abs(qrIK1-qrIK_rtb), abs(qrIK2-qrIK_rtb)));


% Plot robot using Robotics TB
PlotRTB = logical(input('Plot robot animations [1/0]? '));
if ~PlotRTB, return; end
disp(' ');

hf10 = figure(10);
hf10.Color = [1 1 1];
hf10.Position = [200 100 800 800];
set(hf10, 'menubar', 'none');
title('\textbf{SCARA Robot Representation}', ...
      'Interpreter', 'latex','FontSize',14);
xlabel('X','Interpreter', 'latex','FontSize',12);
ylabel('Y','Interpreter', 'latex','FontSize',12);
zlabel('Z','Interpreter', 'latex','FontSize',12);
box on
lgr = legend('Location', 'SouthOutside', 'Orientation', 'horizontal', ...
             'Box', 'off', 'Interpreter', 'latex', 'FontSize', 14);

         
% Plot initial configuration, save fig
lt = 1.75*a1+a2;
pov = [30 36];
title(lgr, 'Initial configuration');
scara.plot(q0, 'workspace', [-lt lt -lt lt 0 d1*1.6], ...
            'notiles', 'noshadow', 'nobase', ...
            'jointdiam', 0.6, 'jointlen', 0.8, ...
            'lightpos', [0.4 0 1], 'fps', 30, ...
            'linkcolor', [1 0.6 0], 'view', pov, ...
            'jointcolor', [0.7 0 1], 'pjointcolor', [0.7 0 1]);
disp(' ');
pause(2);
 
if isempty(dir('Plots')), mkdir('Plots'); end
saveas(hf10,'Plots/FW_Robot_q0.png');


% Plot convergence animation
fpsr = min(300, floor(0.5*length(qout.time)/Tf));
disp('Playing convergence animation...');
title(lgr, 'Playing: Convergence Animation');
scara.plot(qout.data, 'workspace', [-lt lt -lt lt 0 d1*1.6], ...
            'notiles', 'noshadow', 'nobase', ...
            'jointdiam', 0.6, 'jointlen', 0.8, ...
            'lightpos', [0.4 0 1], 'fps', fpsr, ...
            'linkcolor', [1 0.6 0], 'view', pov, ...
            'jointcolor', [0.7 0 1], 'pjointcolor', [0.7 0 1]);
disp('    Done.');
disp(' ');
title(lgr, 'Final configuration');
pause(2);

% Save final configuration plot
saveas(hf10,'Plots/FW_Robot_qr.png');


% Plot trajectory between q0 and qr
disp('Playing normal speed trajectory...');
title(lgr, 'Playing: Trajectory Animation');
qtraj = jtraj([0 0 0 0], qr, linspace(0,1,60));
scara.plot(qtraj, 'workspace', [-lt lt -lt lt 0 d1*1.6], ...
            'notiles', 'noshadow', 'nobase', ...
            'jointdiam', 0.6, 'jointlen', 0.8, ...
            'lightpos', [0.4 0 1], 'fps', 30, ...
            'linkcolor', [1 0.6 0], 'view', pov, ...
            'jointcolor', [0.7 0 1], 'pjointcolor', [0.7 0 1]);
disp('    Done.');
disp(' ');
title(lgr, 'Final configuration');

% Save desired configuration plot
saveas(hf10,'Plots/FW_Robot_qd.png');

if max(abs(qrIK1-qr)) > eStop, qrIKalt = qrIK1;
else, qrIKalt = qrIK2;
end
% Plot alternative final configuration
disp('Plotting alternative final configuration');
title(lgr, 'Final configuration - Alternative');
scara.plot(qrIKalt, 'workspace', [-lt lt -lt lt 0 d1*1.6], ...
            'notiles', 'noshadow', 'nobase', ...
            'jointdiam', 0.6, 'jointlen', 0.8, ...
            'lightpos', [0.4 0 1], 'fps', 30, ...
            'linkcolor', [1 0.6 0], 'view', pov, ...
            'jointcolor', [0.7 0 1], 'pjointcolor', [0.7 0 1]);
disp(' ');

% Save alternative configuration plot
saveas(hf10,'Plots/FW_Robot_qalt.png');
        
disp('Program Ended Correctly')
disp(' ');


