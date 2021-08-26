fig1.Units = 'normalized';
fig1.Position = [0.05 00.05 0.9 0.8];
fig1.Color = [1,1,1];

pov = [-30 20];
All_uSOM = store_state(SOM.coord_ctrl,:);

subplot(2,4,[1,2,5,6])

SOMlength = nxS*nyS;
SOM_ctrl = SOM.coord_ctrl(1:2);
SOM_lowc = S_coord_lc(1:2);
store_pos = store_state(1:3*SOMlength,:);
      
store_x = store_pos(1:SOMlength,:);
limx = [floor(min(store_x(:))*10), ceil(max(store_x(:))*10)]/10;
store_y = store_pos(SOMlength+1:2*SOMlength,:);
limy = [floor(min(store_y(:))*10), ceil(max(store_y(:))*10)]/10;
store_z = store_pos(2*SOMlength+1:3*SOMlength,:);
limz = [floor(min(store_z(:))*10), ceil(max(max(store_z(:)), max(TCP_pos(:,3)))*10)]/10;

wamws = [-0.4 0.4 -0.8 0.2 -0.4 0.8];

plot3(All_uSOM(1,:)',All_uSOM(3,:)',All_uSOM(5,:)','-b');
hold on
plot3(All_uSOM(2,:)',All_uSOM(4,:)',All_uSOM(6,:)','-r');
plot3(store_x(SOM_lowc(1),:)', store_y(SOM_lowc(1),:)', store_z(SOM_lowc(1),:)','color',[1 0.6 0]);
plot3(store_x(SOM_lowc(2),:)', store_y(SOM_lowc(2),:)', store_z(SOM_lowc(2),:)','m');
plot3(Ref_l(:,1),Ref_l(:,2),Ref_l(:,3), '--k');
plot3(Ref_r(:,1),Ref_r(:,2),Ref_r(:,3), '--k');
scatter3(store_x(:,1), store_y(:,1), store_z(:,1), '.b');
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
        fig1.Children(fch).View = [-80 35];
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
plot(time,store_state(S_coord_lc([1 3 5]),:)', 'linewidth',1.5);
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
pa1som = plot(time,store_state(S_coord_lc([2 4 6]),:)', 'linewidth',1.5);
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
               '$x_{SOM}$','$y_{SOM}$', '$z_{SOM}$', 'Ref', ...
               'Orientation','horizontal', 'Interpreter', 'latex');
Lgnd1.Position(1) = 0.75-Lgnd1.Position(3)/2;
Lgnd1.Position(2) = 0.01;
