fig1 = figure(1);
fig1.Units = 'normalized';
fig1.Color = [1,1,1];
fig1.Position = [0.2 0.5 0.5 0.4];

Ts=0.02;

Ref1 = 6;
Ref_l = load(['..\data\trajectories\ref_',num2str(Ref1),'L.csv']);
Ref_r = load(['..\data\trajectories\ref_',num2str(Ref1),'R.csv']);
Ref_lr = [Ref_l, Ref_r]';
Ref_lr = Ref_lr([1,4,2,5,3,6],:);
dR = diff(Ref_lr,1,2);
ddR = diff(dR,1,2);
subplot(11,4,1:4:40)
plot(0:Ts:size(ddR,2)*Ts-Ts, 1e3*ddR([1,3,5],:)','LineWidth',1)
xlim([0 size(ddR,2)*0.02+0.04]);
ylim([-5 3])
grid on
set(gca, 'TickLabelInterpreter', 'latex');
SddR0 = sum(abs(ddR(:)));
MddR0 = 1000*SddR0/size(ddR,2);
xlabel('Time [s]', 'Interpreter', 'latex','FontSize',10)
ylabel('$d^2r/dt^2$ [mm/step$^2$]', 'Interpreter', 'latex','FontSize',10)
title('\textbf{Traj. 0}', 'Interpreter', 'latex','FontSize',10)

Ref1 = 3;
Ref_l = load(['..\data\trajectories\ref_',num2str(Ref1),'L.csv']);
Ref_r = load(['..\data\trajectories\ref_',num2str(Ref1),'R.csv']);
Ref_lr = [Ref_l, Ref_r]';
Ref_lr = Ref_lr([1,4,2,5,3,6],:);
dR = diff(Ref_lr,1,2);
ddR = diff(dR,1,2);
subplot(11,4,2:4:40)
plot(0:Ts:size(ddR,2)*Ts-Ts, 1e3*ddR([1,3,5],:)','LineWidth',1)
xlim([0 size(ddR,2)*0.02+0.04]);
ylim([-5 3])
grid on
set(gca, 'TickLabelInterpreter', 'latex');
SddR1 = sum(abs(ddR(:)));
MddR1 = 1000*SddR1/size(ddR,2);
xlabel('Time [s]', 'Interpreter', 'latex','FontSize',10)
ylabel('$d^2r/dt^2$ [mm/step$^2$]', 'Interpreter', 'latex','FontSize',10)
title('\textbf{Traj. 1}', 'Interpreter', 'latex','FontSize',10)

Ref1 = 13;
NPtT8 = 900/2;
Ref_l = load(['..\data\trajectories\ref_',num2str(Ref1),'L.csv']);
Ref_r = load(['..\data\trajectories\ref_',num2str(Ref1),'R.csv']);
Ref_lr = [Ref_l, Ref_r]';
Ref_lr = Ref_lr([1,4,2,5,3,6],:);
dR = diff(Ref_lr,1,2);
ddR = diff(dR,1,2);
subplot(11,4,3:4:40)
plot(0:Ts:size(ddR,2)*Ts-Ts, 1e3*ddR([1,3,5],:)','LineWidth',1)
xlim([0 size(ddR,2)*0.02+0.02]);
ylim([-5 3])
grid on
set(gca, 'TickLabelInterpreter', 'latex');
SddR2 = sum(sum(abs(ddR(:,1:NPtT8))));
MddR2 = 1000*SddR2/NPtT8;
SddR3 = sum(abs(ddR(:)));
MddR3 = 1000*SddR3/size(ddR,2);
xlabel('Time [s]', 'Interpreter', 'latex','FontSize',10)
ylabel('$d^2r/dt^2$ [mm/step$^2$]', 'Interpreter', 'latex','FontSize',10)
title('\textbf{Trajs. 2+3}', 'Interpreter', 'latex','FontSize',10)
hold on
plot([NPtT8*Ts NPtT8*Ts], ylim,'--k')
hold off

Ref1 = 12;
Ref_l = load(['..\data\trajectories\ref_',num2str(Ref1),'L.csv']);
Ref_r = load(['..\data\trajectories\ref_',num2str(Ref1),'R.csv']);
Ref_lr = [Ref_l, Ref_r]';
Ref_lr = Ref_lr([1,4,2,5,3,6],:);
dR = diff(Ref_lr,1,2);
ddR = diff(dR,1,2);
subplot(11,4,4:4:40)
hpa = plot(0:Ts:size(ddR,2)*Ts-Ts, 1e3*ddR([1,3,5],:)','LineWidth',1);
xlim([0 size(ddR,2)*0.02+0.02]);
ylim([-5 3])
grid on
set(gca, 'TickLabelInterpreter', 'latex');
SddR4 = sum(abs(ddR(:)));
MddR4 = 1000*SddR4/size(ddR,2);
xlabel('Time [s]', 'Interpreter', 'latex','FontSize',10)
ylabel('$d^2r/dt^2$ [mm/step$^2$]', 'Interpreter', 'latex','FontSize',10)
title('\textbf{Traj. 4}', 'Interpreter', 'latex','FontSize',10)

lgnd = legend(hpa, '$X$ \quad', '$Y$ \quad', '$Z$', 'Interpreter', 'latex', ...
    'Orientation', 'horizontal','FontSize',9);
lgnd.Position(1) = 0.5 - lgnd.Position(3)/2;
lgnd.Position(2) = 0.02;
