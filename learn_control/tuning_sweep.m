%% Initialization

clear; clc; close all;

Exps = 5;
filename = 'simlin_traj6';
save_noload = 1;

%Qs=[0:0.05:1,ones(1,20), 0:0.05:1, 0:0.05:0.5, 0:0.1:1];
%Rs=[ones(1,20),1:-0.05:0, 0:0.05:1, 0:0.1:1, 0:0.05:0.5];

Qs=[0:0.05:1,ones(1,20)];
Rs=[ones(1,20),1:-0.05:0];

RMSEs = zeros(length(Qs),1);
TIMEs = zeros(length(Qs),1);


%% Main Loop

if (save_noload==1)
    for wi=1:length(Qs)

        clearvars -except Qs Rs RMSEs TIMEs wi

        W_Q = Qs(wi);
        W_R = Rs(wi);
        simulation_cl_lin;
        RMSEs(wi) = 1000*eRMSEp;
        TIMEs(wi) = tT/nPtRef*1000;

    end
    clearvars -except Qs Rs RMSEs TIMEs filename wi 

    close all;

    save(['RMSEs_',filename,'.mat'], 'RMSEs');
    save(['TIMEs_',filename,'.mat'], 'TIMEs');
else
    load(['RMSEs_',filename,'.mat'], 'RMSEs');
    load(['TIMEs_',filename,'.mat'], 'TIMEs');
end


%% Process

RMSEs = min(RMSEs,10);
fi = RMSEs < min(RMSEs)+0.02*range(RMSEs);
theta = atan2(Rs,Qs);

Qd = Qs./(Rs+Qs);
Rd = 1-Qd;


%% Plots

fig1 = figure(1);
fig1.Color = [1,1,1];
fig1.Units = 'normalized';
fig1.Position = [0.1 0.25 0.25 0.4];

cmj = jet(1000);

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


%scatter(Qd(fi),Rd(fi),[],RMSEs(fi),'filled');


%scatter(theta',RMSEs',[],RMSEs','filled');
%scatter(Qd,RMSEs,[],RMSEs,'filled');
%scatter(Qd(fi),RMSEs(fi),[],RMSEs(fi),'filled');
%scatter(1./ratio(fi),RMSEs(fi),[],RMSEs(fi),'filled');



