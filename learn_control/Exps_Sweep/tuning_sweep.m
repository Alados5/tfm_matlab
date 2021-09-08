%% Initialization

clear; clc; close all;

filename = 'simlin_traj6_noisy';
save_noload = 0;

opts.NTraj = 6;
opts.Ts = 0.02;
opts.Hp = 20;
opts.Wv = 0.15;
opts.nSOM = 4;
opts.nCOM = 4;
opts.nNLM = 10;
opts.ubound = 50*1e-3;
opts.opt_noise = 1;


Qs=[0.05:0.05:1,ones(1,19)];
Rs=[ones(1,19),1:-0.05:0.05];
%Qs=[0:0.05:1,ones(1,20), 0:0.05:1, 0:0.05:0.5, 0:0.1:1];
%Rs=[ones(1,20),1:-0.05:0, 0:0.05:1, 0:0.1:1, 0:0.05:0.5];

RMSEs = zeros(length(Qs),1);
TIMEs = zeros(length(Qs),1);


%% Main Loop

if (save_noload==1)
    for wi=1:length(Qs)

        clearvars -except Qs Rs RMSEs TIMEs filename wi opts

        W_Q = Qs(wi);
        W_R = Rs(wi);
        fprintf(['\nWeights: \t Q=',num2str(W_Q),' \t R=',num2str(W_R),'\n']);
        
        sim_cl_lin_sweep;
        %sim_cl_rtm_sweep;
        RMSEs(wi) = 1000*eRMSEp;
        TIMEs(wi) = tT/nPtRef*1000;

    end
    clearvars -except Qs Rs RMSEs TIMEs filename wi 

    close all;

    save(['data/RMSEs_',filename,'.mat'], 'RMSEs');
    save(['data/TIMEs_',filename,'.mat'], 'TIMEs');
    
else
    load(['data/RMSEs_',filename,'.mat'], 'RMSEs');
    load(['data/TIMEs_',filename,'.mat'], 'TIMEs');
    
end


%% Process

RMSEs = min(RMSEs,10);
fi = RMSEs < min(RMSEs)+0.2*range(RMSEs);

theta = atan2(Rs,Qs);
ratio = Rs./Qs;
Qd = Qs./(Rs+Qs);
Rd = 1-Qd;


%% Plots

cmj = jet(1000);

fig1 = figure(1);
fig1.Color = [1,1,1];
fig1.Units = 'normalized';
fig1.Position = [0.1 0.15 0.45 0.6];

hsp = subplot(2,2,1);
plot([1;1;0],[0;1;1],'--k');
hold on
scatter(Qs,Rs,[],RMSEs,'filled')
hold off
grid on; box on; axis equal;
xlim([0 1.1])
ylim([0 1.1])
xlabel('$Q$','Interpreter','latex','FontSize',10)
ylabel('$R$','Interpreter','latex','FontSize',10)
hsp.FontSize=10;
hsp.TickLabelInterpreter='latex';
colormap(cmj(100:900,:))

hsp = subplot(2,2,2);
plot([1;1;0],[0;1;1],'--k');
hold on
scatter(Qs,Rs,[],TIMEs,'filled')
hold off
grid on; box on; axis equal;
xlim([0 1.1])
ylim([0 1.1])
xlabel('$Q$','Interpreter','latex','FontSize',10)
ylabel('$R$','Interpreter','latex','FontSize',10)
hsp.FontSize=10;
hsp.TickLabelInterpreter='latex';
colormap(cmj(100:900,:))


hsp = subplot(2,2,3);
plot([1;1;0],[0;1;1],'--k');
hold on
scatter(Qs(fi),Rs(fi),[],RMSEs(fi),'filled')
hold off
grid on; box on; axis equal;
xlim([0 1.1])
ylim([0 1.1])
xlabel('$Q$','Interpreter','latex','FontSize',10)
ylabel('$R$','Interpreter','latex','FontSize',10)
hsp.FontSize=10;
hsp.TickLabelInterpreter='latex';
colormap(cmj(100:900,:))

hsp = subplot(2,2,4);
plot([1;1;0],[0;1;1],'--k');
hold on
scatter(Qs(fi),Rs(fi),[],TIMEs(fi),'filled')
hold off
grid on; box on; axis equal;
xlim([0 1.1])
ylim([0 1.1])
xlabel('$Q$','Interpreter','latex','FontSize',10)
ylabel('$R$','Interpreter','latex','FontSize',10)
hsp.FontSize=10;
hsp.TickLabelInterpreter='latex';
colormap(cmj(100:900,:))



fig2 = figure(2);
fig2.Color = [1,1,1];
fig2.Units = 'normalized';
fig2.Position = [0.05 0.15 0.9 0.4];

hsp = subplot(1,3,1);
plot([1;1;0],[0;1;1],'--k');
hold on
scatter(Qs(1:41),Rs(1:41),[],RMSEs(1:41),'filled')
hold off
grid on; box on; axis equal;
xlim([0 1.1])
ylim([0 1.1])
xlabel('$Q$','Interpreter','latex','FontSize',10)
ylabel('$R$','Interpreter','latex','FontSize',10)
hsp.FontSize=11;
hsp.TickLabelInterpreter='latex';
colormap(cmj(100:900,:))
cb1 = colorbar('TickLabelInterpreter','latex');
cb1.Label.String = 'RMSE [mm]';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize = 11;
title('\textbf{All results}', 'Interpreter', 'latex', 'FontSize',14)


hsp = subplot(1,3,2);
plot([1;1;0],[0;1;1],'--k');
hold on
scatter(Qs(fi(1:41)),Rs(fi(1:41)),[],RMSEs(fi(1:41)),'filled')
hold off
grid on; box on; axis equal;
xlim([0 1.1])
ylim([0 1.1])
xlabel('$Q$','Interpreter','latex','FontSize',10)
ylabel('$R$','Interpreter','latex','FontSize',10)
hsp.FontSize=11;
hsp.TickLabelInterpreter='latex';
colormap(cmj(100:900,:))
cb2 = colorbar('TickLabelInterpreter','latex');
cb2.Label.String = 'RMSE [mm]';
cb2.Label.Interpreter = 'latex';
cb2.Label.FontSize = 11;
title('\textbf{Best 2\% of results}', 'Interpreter', 'latex', 'FontSize',14)

hsp = subplot(1,3,3);
plot([1;1;0],[0;1;1],'--k');
hold on
scatter(Qs(fi),Rs(fi),[],RMSEs(fi),'filled')
hold off
grid on; box on; axis equal;
xlim([0 1.1])
ylim([0 1.1])
xlabel('$Q$','Interpreter','latex','FontSize',10)
ylabel('$R$','Interpreter','latex','FontSize',10)
hsp.FontSize=11;
hsp.TickLabelInterpreter='latex';
colormap(cmj(100:900,:))
cb3 = colorbar('TickLabelInterpreter','latex');
cb3.Label.String = 'RMSE [mm]';
cb3.Label.Interpreter = 'latex';
cb3.Label.FontSize = 11;
title('\textbf{Added points with same ratio}', 'Interpreter', 'latex', 'FontSize',14)


%{
fig3 = figure(3);
fig3.Color = [1,1,1];
fig3.Units = 'normalized';
fig3.Position = [0.05 0.15 0.9 0.4];
%}

%scatter(Qd(fi),Rd(fi),[],RMSEs(fi),'filled');


%scatter(theta',RMSEs',[],RMSEs','filled');
%scatter(Qd,RMSEs,[],RMSEs,'filled');
%scatter(Qd(fi),RMSEs(fi),[],RMSEs(fi),'filled');
%scatter(1./ratio(fi),RMSEs(fi),[],RMSEs(fi),'filled');



