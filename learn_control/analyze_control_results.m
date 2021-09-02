clear; clc;

%% Initialization

saveTable = 0;

CtrlData = readtable('LearntCtrl_Data.csv');
WColNames = {'Qx','Qy','Qz','Rx','Ry','Rz'};

% Keep only experiments to analyze
CtrlData = CtrlData(CtrlData.Exps==3, :);
%CtrlData = CtrlData(~sum(CtrlData.Sim==[1,2], 2), :);
CtrlData = CtrlData(~sum(CtrlData.Sim==[0,1], 2), :);

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
