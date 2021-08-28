clear; close all; clc;

% Load all
somstate_00 = load('somstate_traj0_pars0.mat');
somstate_01 = load('somstate_traj0_pars90.mat');
somstate_10 = load('somstate_traj90_pars0.mat');
somstate_11 = load('somstate_traj90_pars90.mat');
somstate_00 = somstate_00.store_state;
somstate_01 = somstate_01.store_state;
somstate_10 = somstate_10.store_state;
somstate_11 = somstate_11.store_state;

somstates = {somstate_00, somstate_01, somstate_10, somstate_11};

nx = 10;
ny = 10;
Nd_S = nx*ny;

coord_nl = [1 nx 1+Nd_S Nd_S+nx 2*Nd_S+1 2*Nd_S+nx];

SOM_lowc = coord_nl(1:2); 
SOM_nd_ctrl = [nx*(ny-1)+1, nx*ny];
SOM_coord_ctrl = [SOM_nd_ctrl, SOM_nd_ctrl+nx*ny, SOM_nd_ctrl+2*nx*ny];

fig3 = figure(3);
fig3.Color = [1,1,1];
fig3.Units = 'normalized';
fig3.Position = [0.2 0.2 0.7 0.5];

for pk = 2:length(somstates)
    store_state = somstates{pk};

    store_pos = store_state(1:3*Nd_S,:);
    store_x = store_pos(1:Nd_S,:);
    limx = [floor(min(store_x(:))*10), ceil(max(store_x(:))*10)]/10;
    store_y = store_pos(Nd_S+1:2*Nd_S,:);
    limy = [floor(min(store_y(:))*10), ceil(max(store_y(:))*10)]/10;
    store_z = store_pos(2*Nd_S+1:3*Nd_S,:);
    limz = [floor(min(store_z(:))*10), ceil(max(store_z(:))*10)]/10;

    All_uSOM = store_state(SOM_coord_ctrl,:);

    subplot(1,3,pk-1);
    plot3(All_uSOM(1:2,:)',All_uSOM(3:4,:)',All_uSOM(5:6,:)');
    hold on
    plot3(store_x(SOM_lowc,:)', store_y(SOM_lowc,:)', store_z(SOM_lowc,:)');
    hold off
    axis equal; box on; grid on;
    xlim(limx);
    ylim(limy);
    zlim(limz);
    set(gca, 'TickLabelInterpreter','latex');
    xlabel('$x$ [m]', 'Interpreter','latex');
    ylabel('$y$ [m]', 'Interpreter','latex');
    zlabel('$z$ [m]', 'Interpreter','latex');
    %title('\textbf{Corner trajectories}','Interpreter','latex');

end









