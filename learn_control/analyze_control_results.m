clear; clc;

%% Initialization

CtrlData = readtable('LearntCtrl_Data.csv');
WColNames = {'Qx','Qy','Qz','Rx','Ry','Rz'};

% Keep only experiments to analyze
CtrlData = CtrlData(CtrlData.Exps==3, :);
CtrlData = CtrlData(~sum(CtrlData.Sim==[1,2], 2), :);
%CtrlData = CtrlData(~sum(CtrlData.Sim==[0,1], 2), :);

% Experiments for analysis: Exps3, Sim0, Traj6, Ts20, Hp25, n4-4, s0
VarColNames = {'Du','Qa','fRw','Wgh'};
KPIColNames = {'RMSE','TOV'};

% Filter unnecessary columns
CtrlData = CtrlData(:, [VarColNames(:)', WColNames(:)', KPIColNames(:)']);
WCols = contains(CtrlData.Properties.VariableNames, WColNames);

AllRwFcns = unique(CtrlData.fRw);
AllWghTypes = unique(CtrlData.Wgh);


%% nCOM & Ts loop

fig1 = figure(1);
fig1.Color = [1,1,1];
fig1.Units = 'normalized';
fig1.Position = [0.0 0.6 0.3 0.3];
fig2 = figure(2);
fig2.Color = [1,1,1];
fig2.Units = 'normalized';
fig2.Position = [0.3 0.6 0.3 0.3];
fig3 = figure(3);
fig3.Color = [1,1,1];
fig3.Units = 'normalized';
fig3.Position = [0.6 0.4 0.4 0.5];

cmj = jet(1000);
cml = lines(4);
spi = 1;

KPIData = table2array(CtrlData(:,{'RMSE','TOV'}));
KPIData01 = (KPIData-min(KPIData))./range(KPIData);

for RFi=AllRwFcns'

    RwFcnData = CtrlData(CtrlData.fRw==RFi,:);


    for WTi=AllWghTypes'

        WghTypeData = RwFcnData(RwFcnData.Wgh==WTi,:);

        X4Data = table2array(WghTypeData(:,{'Du','Qa'}));
        W4Data = table2array(WghTypeData(:,WCols));
        R4Data = table2array(WghTypeData(:,KPIColNames));
        C4Data = KPIData01((CtrlData.Wgh==WTi) & (CtrlData.fRw==RFi), :);

        figure(1);
        hsp = subplot(length(AllRwFcns),length(AllWghTypes),spi);
        scatter(X4Data(:,1),X4Data(:,2),[],C4Data(:,1),'filled');
        box on; axis equal
        xlim([0 1]); ylim([0 1]);
        xticks([0 1]); yticks([0 1])
        xlabel('Using $\Delta u$?','Interpreter','latex','FontSize',10)
        ylabel('Using $Q_a$?','Interpreter','latex','FontSize',10)
        hsp.FontSize=10;
        hsp.TickLabelInterpreter='latex';
        caxis manual
        caxis([0 1]);
        colormap(cmj(100:900,:))

        figure(2);
        hsp = subplot(length(AllRwFcns),length(AllWghTypes),spi);
        scatter(X4Data(:,1),X4Data(:,2),[],C4Data(:,2),'filled');
        box on; axis equal
        xlim([0 1]); ylim([0 1]);
        xticks([0 1]); yticks([0 1])
        xlabel('Using $\Delta u$?','Interpreter','latex','FontSize',10)
        ylabel('Using $Q_a$?','Interpreter','latex','FontSize',10)
        hsp.FontSize=10;
        hsp.TickLabelInterpreter='latex';
        caxis manual
        caxis([0 1]);
        colormap(cmj(100:900,:))

        figure(3);
        subplot(length(AllRwFcns),length(AllWghTypes),spi)
        if(WTi==1)
            plot([1;1;0],[0;1;1],'--k');
            hold on
            hp3=[];
            for r4i=1:4
                hpi=scatter(W4Data(r4i,1), W4Data(r4i,4), [], r4i, 'filled');
                hp3=[hp3;hpi];
            end
            hold off
            axis equal
            xlim([0 1.1]); ylim([0 1.1]);
            grid on; grid minor;
            xlabel('$Q$','Interpreter','latex','FontSize',10)
            ylabel('$R$','Interpreter','latex','FontSize',10)
            set(gca, 'TickLabelInterpreter', 'latex');
            colormap lines;

        elseif(WTi==2)
            scatter3(W4Data(:,1),W4Data(:,2),W4Data(:,3),[],1:4,'filled');
            hold on
            scatter3(W4Data(:,4),W4Data(:,4),W4Data(:,4),[],1:4,'filled');
            plot3([1;1;0;0;1;1],[0;1;1;0;0;0],[1;1;1;1;1;0],'--k');
            plot3([0;0;1;1;0;0],[1;0;0;1;1;1],[0;0;0;0;0;1],'--k');
            plot3([0;0;1;1],[0;0;1;1],[1;0;1;0],'--k');
            hold off
            axis equal
            xlim([0 1.1]); ylim([0 1.1]); zlim([0 1.1]);
            view([110 30]);
            grid on; grid minor;
            xlabel('$x$','Interpreter','latex','FontSize',10)
            ylabel('$y$','Interpreter','latex','FontSize',10)
            zlabel('$z$','Interpreter','latex','FontSize',10)
            set(gca, 'TickLabelInterpreter', 'latex');
            colormap lines;

        elseif(WTi==3)
            scatter3(W4Data(:,1),W4Data(:,2),W4Data(:,3),[],1:4,'filled');
            hold on
            plot3([1;1;0;0;1;1],[0;1;1;0;0;0],[1;1;1;1;1;0],'--k');
            plot3([0;0;1;1;0;0],[1;0;0;1;1;1],[0;0;0;0;0;1],'--k');
            plot3([0;0],[0;0],[1;0],'--k');
            plot3([1;1],[1;1],[1;0],'--k');
            for r4i=1:4
                [~,qrmax]=max([W4Data(r4i,1),-1,-1,W4Data(r4i,4)]);
                plot3([0; W4Data(r4i,qrmax)],[0;W4Data(r4i,qrmax+1)],[0;W4Data(r4i,qrmax+2)],'Color',cml(r4i,:));
            end
            scatter3(W4Data(:,4),W4Data(:,5),W4Data(:,6),[],1:4,'filled');
            hold off
            axis equal
            axlim = max(ceil(max(W4Data(:))*10)/10, 1.1);
            xlim([0 axlim]); ylim([0 axlim]); zlim([0 axlim]);
            view([105 30]);
            grid on; grid minor;
            xlabel('$x$','Interpreter','latex','FontSize',10)
            ylabel('$y$','Interpreter','latex','FontSize',10)
            zlabel('$z$','Interpreter','latex','FontSize',10)
            set(gca, 'TickLabelInterpreter', 'latex');
            colormap lines;

        end

        if(spi==1)
            lgnd = legend(hp3, '$u,Q_k$\quad','$\Delta u,Q_k$\quad', ...
                    '$u,Q_a$\quad','$\Delta u,Q_a$\quad', ...
                    'orientation','horizontal','Interpreter','latex');
            lgnd.Position(1) = 0.5-lgnd.Position(3)/2;
            lgnd.Position(2) = 0.01;
        end

        spi = spi+1;

    end
end


%% Add colorbars

NTicks = 6;

figure(1);
for spi=1:length(AllRwFcns)*length(AllWghTypes)

    hsp = subplot(length(AllRwFcns),length(AllWghTypes),spi);
    hsp.Position(1)=hsp.Position(1)-0.05;

end
cb1ticks = min(KPIData(:,1))+linspace(0,1,NTicks)*range(KPIData(:,1));
cb1ticks = round(cb1ticks*100)/100;
cb1 = colorbar('Position',[0.89 0.149 0.02 0.776], ...
               'TickLabels',cb1ticks, 'TickLabelInterpreter','latex');
cb1.Label.String = 'RMSE [mm]';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize = 10;

figure(2);
for spi=1:length(AllRwFcns)*length(AllWghTypes)

    hsp = subplot(length(AllRwFcns),length(AllWghTypes),spi);
    hsp.Position(1)=hsp.Position(1)-0.05;

end
cb2ticks = min(KPIData(:,2))+linspace(0,1,NTicks)*range(KPIData(:,2));
cb2ticks = round(cb2ticks*100)/100;
cb2 = colorbar('Position',[0.89 0.149 0.02 0.776], ...
               'TickLabels',cb2ticks, 'TickLabelInterpreter','latex');
cb2.Label.String = '$t_c/T_s$';
cb2.Label.Interpreter = 'latex';
cb2.Label.FontSize = 10;


%% For Exps4

clear; clc;

loadtraj = 'traj3_ts20_hp25_ns4_nc4';

load(['Exps4/0_LIN_Det/',loadtraj,'/TH_0-2.mat'])
load(['Exps4/0_LIN_Det/',loadtraj,'/RW_0-2.mat'])

cmj = jet(1000);

Data0 = [TH(:,:,1)', RW(:,:,1)'];
Data1 = [TH(:,:,2)', RW(:,:,2)'];

Data0f = Data0(Data0(:,3)~=-10, :);
Data0f(:,3) = (Data0f(:,3)-min(Data0f(:,3)))/range(Data0f(:,3));
Data0f = Data0f(Data0f(:,3)>0.98, :);
Data0f(:,3) = (Data0f(:,3)-min(Data0f(:,3)))/range(Data0f(:,3));

Data1f = Data1(Data1(:,3)~=-10, :);
Data1f_abs = abs(Data1f(:,3));
Data1f(:,3) = (Data1f(:,3)-min(Data1f(:,3)))/range(Data1f(:,3));
threshold = 0.5;
Data1f_abs = Data1f_abs(Data1f(:,3)>threshold);
Data1f = Data1f(Data1f(:,3)>threshold, :);
Data1f(:,3) = (Data1f(:,3)-min(Data1f(:,3)))/range(Data1f(:,3));

theta = atan2(Data1f(:,2), Data1f(:,1));
[theta, id] = sort(theta);
data1t = [theta Data1f_abs(id)];

ratio = Data1f(:,2)./Data1f(:,1);
[ratio, id] = sort(ratio);
data1o = [ratio Data1f_abs(id)];

[Rsort, id] = sort(Data1f(:,2));
data1r = [Rsort Data1f_abs(id)];
data1r(data1r(:,1)==1,:)=[];

data1n = [data1r(:,1), data1r(:,2) - min(data1r(:,2))];


% Plot rainbow scale for one experiment
fig4 = figure(4);
fig4.Color = [1,1,1];
fig4.Units = 'normalized';
fig4.Position = [0.05 0.1 0.9 0.4];

hsp=subplot(1,3,1);
plot([1;1;0],[0;1;1],'--k');
hold on
scatter(Data0(:,1), Data0(:,2), [], -Data0(:,3), 'filled')
hold off
axis equal; grid on
xlim([0 1.1])
ylim([0 1.1])
xlabel('$Q$','Interpreter','latex','FontSize',10)
ylabel('$R$','Interpreter','latex','FontSize',10)
hsp.FontSize=11;
hsp.TickLabelInterpreter='latex';
colormap(cmj(100:900,:))
cb7 = colorbar('TickLabelInterpreter','latex');
cb7.Label.String = 'RMSE [mm]';
cb7.Label.Interpreter = 'latex';
cb7.Label.FontSize = 11;
title('\textbf{First epoch results}', 'Interpreter', 'latex', 'FontSize',14)

hsp=subplot(1,3,2);
plot([1;1;0],[0;1;1],'--k');
hold on
scatter(Data1(:,1), Data1(:,2), [], -Data1(:,3), 'filled')
hold off
axis equal; grid on
xlim([0 1.1])
ylim([0 1.1])
xlabel('$Q$','Interpreter','latex','FontSize',10)
ylabel('$R$','Interpreter','latex','FontSize',10)
hsp.FontSize=11;
hsp.TickLabelInterpreter='latex';
colormap(cmj(100:900,:))
cb7 = colorbar('TickLabelInterpreter','latex');
cb7.Label.String = 'RMSE [mm]';
cb7.Label.Interpreter = 'latex';
cb7.Label.FontSize = 11;
title('\textbf{Second epoch results}', 'Interpreter', 'latex', 'FontSize',14)

hsp=subplot(1,3,3);
plot([1;1;0],[0;1;1],'--k');
hold on
%scatter(Data0f(:,1), Data0f(:,2), [], 1-Data0f(:,3), 'filled')
scatter(Data1f(:,1), Data1f(:,2), [], Data1f_abs, 'filled')
hold off
axis equal; grid on
xlim([0 1.1])
ylim([0 1.1])
xlabel('$Q$','Interpreter','latex','FontSize',10)
ylabel('$R$','Interpreter','latex','FontSize',10)
hsp.FontSize=11;
hsp.TickLabelInterpreter='latex';
colormap(cmj(100:900,:))
cb7 = colorbar('TickLabelInterpreter','latex');
cb7.Label.String = 'RMSE [mm]';
cb7.Label.Interpreter = 'latex';
cb7.Label.FontSize = 11;
title('\textbf{Best 2\% of second epoch results}', 'Interpreter', 'latex', 'FontSize',14)



% Compare experiments
fig5 = figure(5);
fig5.Color = [1,1,1];
fig5.Units = 'normalized';
fig5.Position = [0.05 0.3 0.7 0.4];

Trajs = [6,3,13,12,16];
NFiles = length(Trajs);

hsp=subplot(1,1,1);
for f=1:NFiles
    filename = ['traj',num2str(Trajs(f)),'_ts20_hp25_ns4_nc4'];
    
    load(['Exps4/0_LIN_Det/',filename,'/TH_0-2.mat'])
    load(['Exps4/0_LIN_Det/',filename,'/RW_0-2.mat'])
    
    Data1 = [TH(:,:,2)', RW(:,:,2)'];

    Data1f = Data1(Data1(:,3)~=-10, :);
    Data1f_abs = abs(Data1f(:,3));
    Data1f(:,3) = (Data1f(:,3)-min(Data1f(:,3)))/range(Data1f(:,3));
    threshold = 0.5;
    Data1f_abs = Data1f_abs(Data1f(:,3)>threshold);
    Data1f = Data1f(Data1f(:,3)>threshold, :);
    Data1f(:,3) = (Data1f(:,3)-min(Data1f(:,3)))/range(Data1f(:,3));

    [Rsort, id] = sort(Data1f(:,2));
    data1r = [Rsort Data1f_abs(id)];
    data1r(data1r(:,1)==1,:)=[];

    data1n = [data1r(:,1), (data1r(:,2) - min(data1r(:,2)))/range(data1r(:,2))];
    
    if (f>1), hold on; end
    plot(data1n(:,1), data1n(:,2), 'LineWidth',1);

end

hold off
grid on
xlim([0 1])
xlabel('Ratio $R/Q$','Interpreter','latex','FontSize',10)
ylabel('Normalized RMSE','Interpreter','latex','FontSize',10)
hsp.FontSize=11;
hsp.TickLabelInterpreter='latex';
title('\textbf{Results for all the analyzed trajectories}', 'Interpreter', 'latex', 'FontSize',14)

legend('Traj. 0', 'Traj. 1', 'Traj. 2', 'Traj. 3', 'Traj. 4', ...
       'Location','SouthEast','Interpreter','latex');
   
   
%% For Exps5

clear; clc;

TuningTss = [10 15 20 25];
TuningHps = [15 20 25 30];
NFiles = length(TuningTss)*length(TuningHps);

AllMW = zeros(NFiles,2);
AllRW = zeros(NFiles,1);
AllBD = zeros(NFiles,2);
AllHZ = zeros(NFiles,2);

%{
fig6 = figure(6);
fig6.Color = [1,1,1];
fig6.Units = 'normalized';
fig6.Position = [0.05 0.3 0.7 0.4];
%}

cmj = jet(1000);

hsp=subplot(1,1,1);
spi=1;
for Tsi=1:length(TuningTss)
    for Hpi=1:length(TuningHps)
    
        filename = ['traj6_ts',num2str(TuningTss(Tsi)), ...
                    '_hp',num2str(TuningHps(Hpi)),'_ns4_nc4'];

        load(['Exps5/0_LIN_Tuning_Det/',filename,'/MW_0-5.mat'])
        load(['Exps5/0_LIN_Tuning_Det/',filename,'/RWMW_0-5.mat'])
        load(['Exps5/0_LIN_Tuning_Det/',filename,'/BD_0-5.mat'])
        
        MW = MW(:,:,end)';
        RWMW = RWMW(end);
        BD = BD(:,:,end);
        
        AllHZ(spi,:) = [TuningTss(Tsi) TuningHps(Hpi)];
        AllMW(spi,:) = MW;
        AllRW(spi,:) = -RWMW;
        AllBD(spi,:) = BD;
        
        spi=spi+1;

    end

end











