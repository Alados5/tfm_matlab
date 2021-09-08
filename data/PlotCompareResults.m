% Plot 3D and 2D evolutions in a single figure
%   given a store_state .mat file to load

clear; clc;

filename1 = 'somstate_rtm_noise1.mat';
filename2 = 'somstate_rtm_noise2.mat';
Ref1 = 13;
Ref2 = 13;
Ts = 0.02; %0.01;
nSOM = 4;


somstate1 = load(['results\',filename1]);
somstate1 = somstate1.store_somstate;
somstate2 = load(['results\',filename2]);
somstate2 = somstate2.store_somstate;

nPt1 = size(somstate1,2);
nPt2 = size(somstate2,2);
time1 = 0:Ts:nPt1*Ts-Ts;
time2 = 0:Ts:nPt2*Ts-Ts;

Ref1_l = load(['trajectories\ref_',num2str(Ref1),'L.csv']);
Ref1_r = load(['trajectories\ref_',num2str(Ref1),'R.csv']);
Ref2_l = load(['trajectories\ref_',num2str(Ref2),'L.csv']);
Ref2_r = load(['trajectories\ref_',num2str(Ref2),'R.csv']);
Ref1_lr = [Ref1_l, Ref1_r]';
Ref1_lr = Ref1_lr([1,4,2,5,3,6],:);
Ref2_lr = [Ref2_l, Ref2_r]';
Ref2_lr = Ref2_lr([1,4,2,5,3,6],:);

Nd_S = nSOM*nSOM;
coord_nl = [1 nSOM 1+Nd_S Nd_S+nSOM 2*Nd_S+1 2*Nd_S+nSOM];
nd_ctrl = [nSOM*(nSOM-1)+1, nSOM*nSOM];
coord_ctrl = [nd_ctrl, nd_ctrl+nSOM*nSOM, nd_ctrl+2*nSOM*nSOM];

u_SOM1 = somstate1(coord_ctrl,:);
u_SOM2 = somstate2(coord_ctrl,:);
u_lin1 = diff(u_SOM1,1,2);
u_lin2 = diff(u_SOM2,1,2);
u_lin1(5:6,1:20) = 0;
u_lin2(5:6,1:20) = 0;
Du_1 = diff(u_lin1,1,2);
Du_2 = diff(u_lin2,1,2);

Dulim = [min([Du_1(:); Du_2(:)]), max([Du_1(:); Du_2(:)])];
Dulim = round(10000*Dulim)/10;

err1 = somstate1(coord_nl,:) - Ref1_lr;
err1(5:6,1:20) = 0;
err2 = somstate2(coord_nl,:) - Ref2_lr;
err2(5:6,1:20) = 0;

dR1 = diff(Ref1_lr,1,2);
ddR1 = diff(dR1,1,2);
dR2 = diff(Ref2_lr,1,2);
ddR2 = diff(dR2,1,2);

store_pos1 = somstate1(1:3*Nd_S,:);
store_pos2 = somstate2(1:3*Nd_S,:);

store_x1 = store_pos1(1:Nd_S,:);
store_y1 = store_pos1(Nd_S+1:2*Nd_S,:);
store_z1 = store_pos1(2*Nd_S+1:3*Nd_S,:);
store_x2 = store_pos2(1:Nd_S,:);
store_y2 = store_pos2(Nd_S+1:2*Nd_S,:);
store_z2 = store_pos2(2*Nd_S+1:3*Nd_S,:);
limx1 = [floor(min(store_x1(:))*10), ceil(max(store_x1(:))*10)]/10;
limy1 = [floor(min(store_y1(:))*10), ceil(max(store_y1(:))*10)]/10;
limz1 = [floor(min(store_z1(:))*10), ceil(max(store_z1(:))*10)]/10;
limx2 = [floor(min(store_x2(:))*10), ceil(max(store_x2(:))*10)]/10;
limy2 = [floor(min(store_y2(:))*10), ceil(max(store_y2(:))*10)]/10;
limz2 = [floor(min(store_z2(:))*10), ceil(max(store_z2(:))*10)]/10;
limx = [min(limx1(1), limx2(1)), max(limx1(2), limx2(2))];
limy = [min(limy1(1), limy2(1)), max(limy1(2), limy2(2))];
limz = [min(limz1(1), limz2(1)), max(limz1(2), limz2(2))];

pov = [-60 25]; %[-30 20];

fig1 = figure(1);
fig1.Units = 'normalized';
fig1.Position = [0.05 0.05 0.9 0.7];
fig1.Color = [1,1,1];

subplot(2,4,[1,2,5,6])

plot3(u_SOM1(1,:)',u_SOM1(3,:)',u_SOM1(5,:)','-b');
hold on
plot3(u_SOM1(2,:)',u_SOM1(4,:)',u_SOM1(6,:)','-r');
plot3(store_x1(coord_nl(1),:)', store_y1(coord_nl(1),:)', store_z1(coord_nl(1),:)','color',[1 0.6 0]);
plot3(store_x1(coord_nl(2),:)', store_y1(coord_nl(2),:)', store_z1(coord_nl(2),:)','-m');

plot3(u_SOM2(1,:)',u_SOM2(3,:)',u_SOM2(5,:)','--b');
plot3(u_SOM2(2,:)',u_SOM2(4,:)',u_SOM2(6,:)','--r');
plot3(store_x2(coord_nl(1),:)', store_y2(coord_nl(1),:)', store_z2(coord_nl(1),:)','--','color',[1 0.6 0]);
plot3(store_x2(coord_nl(2),:)', store_y2(coord_nl(2),:)', store_z2(coord_nl(2),:)','--m');
hold off
axis equal; box on; grid on;
xlim(limx);
ylim(limy);
zlim(limz);
set(gca, 'TickLabelInterpreter','latex');
xlabel('$X$ [m]', 'Interpreter','latex');
ylabel('$Y$ [m]', 'Interpreter','latex');
zlabel('$Z$ [m]', 'Interpreter','latex');
title('\textbf{3D Evolution of all corners}', 'Interpreter', 'latex')


for fch=1:length(fig1.Children)
    if isa(fig1.Children(fch),'matlab.graphics.axis.Axes')
        fig1.Children(fch).View = pov;
    end
end

%{
subplot(2,4,3)
plot(time1(1:end-2), 1e3*Du_1([1,3,5],:)', 'linewidth',1);
title('\textbf{Slew rate evolution (1)}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('$\Delta u$ [mm]', 'Interpreter', 'latex')
xlim([0 time1(end)])
ylim(Dulim);
set(gca, 'TickLabelInterpreter', 'latex');

subplot(2,4,4)
plot(time2(1:end-2), 1e3*Du_2([1,3,5],:)', 'linewidth',1);
title('\textbf{Slew rate evolution (2)}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('$\Delta u$ [mm]', 'Interpreter', 'latex')
xlim([0 time1(end)])
ylim(Dulim);
set(gca, 'TickLabelInterpreter', 'latex');
%}

%{x
subplot(2,4,3)
plot(time1,somstate1(coord_ctrl([1 3 5]),:)', 'linewidth',1.5);
title('\textbf{Upper left corner evolution (1)}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time1(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(2,4,4)
plot(time1,somstate2(coord_ctrl([1 3 5]),:)', 'linewidth',1.5);
title('\textbf{Upper left corner evolution (2)}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time1(end)])
set(gca, 'TickLabelInterpreter', 'latex');
%}


subplot(2,4,7)
plot(time1,somstate1(coord_nl([1 3 5]),:)', 'linewidth',1.5);
hold on
plot(time1,Ref1_l, ':k', 'linewidth',1.2);
hold off
title('\textbf{Lower left corner evolution (1)}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time1(end)])
set(gca, 'TickLabelInterpreter', 'latex');

subplot(2,4,8)
pa1som = plot(time2,somstate2(coord_nl([1 3 5]),:)', 'linewidth',1.5);
hold on
pa1ref = plot(time2,Ref2_l, ':k', 'linewidth',1.2);
hold off
title('\textbf{Lower left corner evolution (2)}', 'Interpreter', 'latex')
grid on
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Position [m]', 'Interpreter', 'latex')
xlim([0 time1(end)])
set(gca, 'TickLabelInterpreter', 'latex');

Lgnd1 = legend([pa1som', pa1ref(1)], ...
               '$X$','$Y$', '$Z$', '$r$', ...
               'Orientation','horizontal', 'Interpreter', 'latex');
Lgnd1.Position(1) = 0.75-Lgnd1.Position(3)/2;
Lgnd1.Position(2) = 0.01;

%{
title('\textbf{Lower left corner evolution ({\boldmath$\bar{u}=5$} mm)}', 'Interpreter', 'latex', 'FontSize',12)
title('\textbf{Lower left corner evolution ({\boldmath$\bar{u}=50$} mm)}', 'Interpreter', 'latex', 'FontSize',12)
title('\textbf{Upper left corner evolution ({\boldmath$\bar{u}=5$} mm)}', 'Interpreter', 'latex', 'FontSize',12)
title('\textbf{Upper left corner evolution ({\boldmath$\bar{u}=50$} mm)}', 'Interpreter', 'latex', 'FontSize',12)
%}
