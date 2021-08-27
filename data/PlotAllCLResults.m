% Plot 3D and 2D evolutions in a single figure
%   given a store_state .mat file to load

clear; clc;

filename = 'somstate_rotproof.mat';
NTraj = 16;
Ts = 0.01;


load(['results\',filename]);
Ref_l = load(['trajectories\ref_',num2str(NTraj),'L.csv']);
Ref_r = load(['trajectories\ref_',num2str(NTraj),'R.csv']);
nPtRef = size(Ref_l,1);
time = 0:Ts:nPtRef*Ts-Ts;

nx = 10;
ny = 10;
Nd_S = nx*ny;

coord_nl = [1 nx 1+Nd_S Nd_S+nx 2*Nd_S+1 2*Nd_S+nx];

SOM_nd_ctrl = [nx*(ny-1)+1, nx*ny];
SOM_coord_ctrl = [SOM_nd_ctrl, SOM_nd_ctrl+nx*ny, SOM_nd_ctrl+2*nx*ny];

pov = [-30 20];
All_uSOM = store_state(SOM_coord_ctrl,:);

SOMlength = nx*ny;
SOM_ctrl = SOM_coord_ctrl(1:2);
SOM_lowc = coord_nl(1:2); 
store_pos = store_state(1:3*SOMlength,:);
      
store_x = store_pos(1:SOMlength,:);
limx = [floor(min(store_x(:))*10), ceil(max(store_x(:))*10)]/10;
store_y = store_pos(SOMlength+1:2*SOMlength,:);
limy = [floor(min(store_y(:))*10), ceil(max(store_y(:))*10)]/10;
store_z = store_pos(2*SOMlength+1:3*SOMlength,:);
limz = [floor(min(store_z(:))*10), ceil(max(store_z(:))*10)]/10;

fig1 = figure(1);
fig1.Units = 'normalized';
fig1.Position = [0.05 0.05 0.9 0.7];
fig1.Color = [1,1,1];

subplot(2,4,[1,2,5,6])

plot3(All_uSOM(1,:)',All_uSOM(3,:)',All_uSOM(5,:)','-b');
hold on
plot3(All_uSOM(2,:)',All_uSOM(4,:)',All_uSOM(6,:)','-r');
plot3(store_x(SOM_lowc(1),:)', store_y(SOM_lowc(1),:)', store_z(SOM_lowc(1),:)','color',[1 0.6 0]);
plot3(store_x(SOM_lowc(2),:)', store_y(SOM_lowc(2),:)', store_z(SOM_lowc(2),:)','m');
plot3(Ref_l(:,1),Ref_l(:,2),Ref_l(:,3), '--k');
plot3(Ref_r(:,1),Ref_r(:,2),Ref_r(:,3), '--k');
scatter3(store_x(:,1), store_y(:,1), store_z(:,1), [],[0.6 0.6 0.6], '.');
hold off
axis equal; box on; grid on;
xlim(limx);
ylim(limy);
zlim(limz);
set(gca, 'TickLabelInterpreter','latex');
xlabel('$x$ [m]', 'Interpreter','latex');
ylabel('$y$ [m]', 'Interpreter','latex');
zlabel('$z$ [m]', 'Interpreter','latex');
title('\textbf{3D Evolution of all corners}', 'Interpreter', 'latex')


for fch=1:length(fig1.Children)
    if isa(fig1.Children(fch),'matlab.graphics.axis.Axes')
        fig1.Children(fch).View = [-81 27];
    end
end

subplot(2,4,3)
plot(time, All_uSOM([1,3,5],:)','linewidth',1.5)
title('\textbf{Left upper corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(2,4,4)
plot(time, All_uSOM([2,4,6],:)','linewidth',1.5);
title('\textbf{Right upper corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(2,4,7)
plot(time,store_state(coord_nl([1 3 5]),:)', 'linewidth',1.5);
hold on
plot(time,Ref_l, '--k', 'linewidth',1.2);
hold off
title('\textbf{Left lower corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(2,4,8)
pa1som = plot(time,store_state(coord_nl([2 4 6]),:)', 'linewidth',1.5);
hold on
pa1ref = plot(time,Ref_r, '--k', 'linewidth',1.2);
hold off
title('\textbf{Right lower corner}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time(end)])
set(gca, 'TickLabelInterpreter', 'latex');

Lgnd1 = legend([pa1som' pa1ref(1)], ...
               '$x_{SOM}$','$y_{SOM}$', '$z_{SOM}$', '$r$', ...
               'Orientation','horizontal', 'Interpreter', 'latex');
Lgnd1.Position(1) = 0.75-Lgnd1.Position(3)/2;
Lgnd1.Position(2) = 0.01;
