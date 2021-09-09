clear; clc; close all;

addpath('../required_files/cloth_model_New_L')

% Cloth mesh
pos = create_lin_mesh(0.3, 10, [0 -0.4 0.1], 0);

% Cloth upper corners
u_SOM = [-0.1500 0.1500 -0.4000 -0.4000 0.2500 0.2500]';

BaseL = 0.3;

% Get Cloth orientation (rotation matrix)
cloth_x = u_SOM([2 4 6]) - u_SOM([1 3 5]);
cloth_y = [-cloth_x(2) cloth_x(1) 0]';
cloth_z = cross(cloth_x,cloth_y);

cloth_x = cloth_x/norm(cloth_x);
cloth_y = cloth_y/norm(cloth_y);
cloth_z = cloth_z/norm(cloth_z);
Rcloth = [cloth_x cloth_y cloth_z];
Rtcp = [cloth_y cloth_x -cloth_z];

TCP_pos = (u_SOM([1,3,5]) + u_SOM([2,4,6]))/2 + [0;0;0.09];
TCP_rot = Rtcp;
TCP_Tm  = [TCP_rot TCP_pos;
          [ 0 0 0       1]];
EndLT = 0.35*TCP_rot*BaseL+TCP_pos;

      
% Robot-Camera transform   
TRC = [1.0,  0.0,  0.0,  0.1;
       0.0,  0.0,  1.0, -1.4;
       0.0, -1.0,  0.0,  0.2;
       0.0,  0.0,  0.0,  1.0];
RRC = TRC(1:3,1:3);
pRC = TRC(1:3,4);
EndBC = RRC*BaseL+pRC;
EndLC = 0.8*RRC*BaseL+pRC;

% Initialize WAM model
run("../with_robot/init_WAM.m");
q0 = wam.ikine(TCP_Tm, 'q0',qref);


wamws = [-0.6 0.6 -1.6 0.4 -0.4 0.8];
wampov = [-52 16];
WAMbaseC = 0.8*[0.8 0.4 0];
WAMbaseA = 0.5;



fig1 = figure(1);
fig1.Color = [1,1,1];
fig1.Units = 'normalized';
fig1.Position = [0.05 0.1 0.35 0.50];

% Plot Base
hold on
fill3([-0.1 0.1 0.1 -0.1],[-0.1 -0.1 0.1 0.1],[0 0 0 0],WAMbaseC,'FaceAlpha',WAMbaseA)
fill3([-0.1 -0.1 -0.1 -0.1],[-0.1 -0.1 0.1 0.1],[-0.5 0 0 -0.5],WAMbaseC,'FaceAlpha',WAMbaseA)
fill3([-0.1 -0.1 0.1 0.1],[-0.1 -0.1 -0.1 -0.1],[-0.5 0 0 -0.5],WAMbaseC,'FaceAlpha',WAMbaseA)
fill3([0.1 0.1 0.1 0.1],[-0.1 -0.1 0.1 0.1],[-0.5 0 0 -0.5],WAMbaseC,'FaceAlpha',WAMbaseA)
fill3([-0.1 -0.1 0.1 0.1],[0.1 0.1 0.1 0.1],[-0.5 0 0 -0.5],WAMbaseC,'FaceAlpha',WAMbaseA)
hold off

% Plot WAM
wam.plot(q0, 'workspace', wamws, ...
                         'notiles', 'noshadow', 'nobase', ...
                         'jointdiam', 0.6, 'jointlen', 0.8, ...
                         'lightpos', [-0.5 -0.5 1], 'fps', 30, ...
                         'linkcolor', [1 0.6 0], 'view', wampov, ...
                         'jointcolor', [0.7 0 1], 'pjointcolor', [0.7 0 1]);

%
hold on
plot3([0 0.8*BaseL],[0 0],[0 0],'r','LineWidth',2);
plot3([0 0],[0 0.8*BaseL],[0 0],'g','LineWidth',2);
plot3([0 0],[0 0],[0 0.8*BaseL],'b','LineWidth',2);
arrow3([0,0,0],[BaseL,0,0],'r')
arrow3([0,0,0],[0,BaseL,0],'g')
arrow3([0,0,0],[0,0,BaseL],'b')

scatter3(pRC(1),pRC(2),pRC(3),200,'ok','filled');
plot3([pRC(1), pRC(1)], [pRC(2), pRC(2)-0.2], ...
      [pRC(3), -0.4],'k','LineWidth',1)
plot3([pRC(1), pRC(1)+0.3*sqrt(3)/3], [pRC(2), pRC(2)+0.1], ...
      [pRC(3), -0.4],'k','LineWidth',1)
plot3([pRC(1), pRC(1)-0.3*sqrt(3)/3], [pRC(2), pRC(2)+0.1], ...
      [pRC(3), -0.4],'k','LineWidth',1)


plot3([pRC(1), EndLC(1,1)], [pRC(2), EndLC(2,1)], ...
      [pRC(3), EndLC(3,1)],'r','LineWidth',2)
plot3([pRC(1), EndLC(1,2)], [pRC(2), EndLC(2,2)], ...
      [pRC(3), EndLC(3,2)],'g','LineWidth',2)
plot3([pRC(1), EndLC(1,3)], [pRC(2), EndLC(2,3)], ...
      [pRC(3), EndLC(3,3)],'b','LineWidth',2)
arrow3(pRC',EndBC(:,1)','r')
arrow3(pRC',EndBC(:,2)','g')
arrow3(pRC',EndBC(:,3)','b')

scatter3(pos(:,1),pos(:,2),pos(:,3),10,0.5*[1 1 1],'filled')
plot3([TCP_pos(1), EndLT(1,1)], [TCP_pos(2), EndLT(2,1)], ...
      [TCP_pos(3), EndLT(3,1)],'r','LineWidth',2)
plot3([TCP_pos(1), EndLT(1,2)], [TCP_pos(2), EndLT(2,2)], ...
      [TCP_pos(3), EndLT(3,2)],'g','LineWidth',2)
plot3([TCP_pos(1), EndLT(1,3)], [TCP_pos(2), EndLT(2,3)], ...
      [TCP_pos(3), EndLT(3,3)],'b','LineWidth',2)

hold off

set(gca, 'FontSize', 11);
set(gca, 'TickLabelInterpreter','latex');
xlabel('$X$ [m]', 'Interpreter','latex');
ylabel('$Y$ [m]', 'Interpreter','latex');
zlabel('$Z$ [m]', 'Interpreter','latex');


