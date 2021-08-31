clear; clc;

%% Initialization

saveTable = 0;

MdlData = readtable('LearntMdl_Data.csv');
ThCols = contains(MdlData.Properties.VariableNames, 'Th_');
ThColNames = MdlData.Properties.VariableNames(ThCols);
MainColNames = {'ExpSetN','NExp','NTrial','Ts','nSOM','nCOM'};

% Filter unnecessary columns
MdlData = MdlData(:, [MainColNames(:)', ThColNames(:)', {'Rwd'}]);
ThCols = contains(MdlData.Properties.VariableNames, 'Th_');

% Remove exps and trajectories not fit for learning
MdlData = MdlData(~sum(MdlData.ExpSetN==[0:2,11], 2), :);
MdlData = MdlData(~sum(MdlData.NExp==[1:7,10,12:13], 2), :);
%MdlData = MdlData(~sum(MdlData.NTrial==3, 2), :);
MdlData(MdlData.Rwd<-1200,'Rwd')={-1200};

AllExpSets = unique(MdlData.ExpSetN);
AllCOMSizes = unique(MdlData.nCOM);
AllTrajs = unique(MdlData.NExp);
AllTs = unique(MdlData.Ts);


%% ExpSetN & Ts loop

%{
fig1 = figure(1);
fig1.Color = [1,1,1];
fig1.Units = 'normalized';
fig2 = figure(2);
fig2.Color = [1,1,1];
fig2.Units = 'normalized';
fig2.Position(3:4) = [0.4 0.35];

TableMeansE = [];
TableMediansE = [];
TableWavgsE = [];

colormap jet;
TsColors = [[0 0.6 0]; [0.1 0.4 1]; [0.8 0 1]; [1 0 0]];

for ESi=AllExpSets'
    
    ExpSetData = MdlData(MdlData.ExpSetN==ESi,:);
    
    for TSi=AllTs'
        
        TsMdlData = ExpSetData(ExpSetData.Ts==TSi,:);
        ThTsMdlData = table2array(TsMdlData(:,ThCols))';
        RwTsMdlData = 1 - (TsMdlData.Rwd - min(TsMdlData.Rwd))/(range(TsMdlData.Rwd)+eps);
        PlotX = kron(ones(size(TsMdlData,1),1), (1:7)');
        PlotC = kron(RwTsMdlData, ones(7,1));
        
        if size(ThTsMdlData,2)>0
            dw = REPSupdate(TsMdlData.Rwd);
            
            ThTsMdlMean = mean(ThTsMdlData,2);
            ThTsMdlMedian = median(ThTsMdlData,2);
            ThTsMdlWavg = sum(ThTsMdlData.*dw,2)/sum(dw);
            
            TableMeansE = [TableMeansE; ESi TSi ThTsMdlMean'];
            TableMediansE = [TableMediansE; ESi TSi ThTsMdlMedian'];
            TableWavgsE = [TableWavgsE; ESi TSi ThTsMdlWavg'];
                        
            figure(2);
            subplot(1,9,1:3)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            %scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            scatter((1:7)', ThTsMdlWavg, 50, 'om', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([0.5,3.5]);
            xticklabels({'$k_x$','$k_y$','$k_z$'});
            ylabel('Value [N/m]','Interpreter','latex','FontSize',10);
            subplot(1,9,5:7)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            %scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            scatter((1:7)', ThTsMdlWavg, 50, 'om', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([3.5,6.5]);
            xticklabels({'$b_x$','$b_y$','$b_z$'});
            ylabel('Value [Ns/m]','Interpreter','latex','FontSize',10);
            subplot(1,9,9)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            %scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            scatter((1:7)', ThTsMdlWavg, 50, 'om', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([6.5,7.5]);
            xticklabels({'','$\Delta l_{0z}$',''});
            ylabel('Value [m]','Interpreter','latex','FontSize',10);
            
            %{
            figure(1);
            subplot(length(AllTs), length(AllExpSets), spi)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            if size(ThTsMdlData,2)>1
                hold on
                boxplot(ThTsMdlData');
                hold off
            end
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            %}
            %{
            figure(1);
            hold on
            scatter3(ESi*ones(size(TsMdlData,1),1),TSi*ones(size(TsMdlData,1),1),TsMdlData.Rwd,25,RwTsMdlData,'filled');
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            zlim([-1000 0]);
            %}
            figure(1);
            subplot(1,2,1)
            hold on
            scatter(ESi*ones(size(TsMdlData,1),1), TsMdlData.Rwd,25,TsColors(AllTs==TSi,:).*ones(size(TsMdlData,1),1),'filled');
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            ylim([min(max(-1000, min(TsMdlData.Rwd)), min(ylim)) 0]);
            xlim([min(AllExpSets)-1, max(AllExpSets)+1]);
            subplot(1,2,2)
            hold on
            scatter(TSi*ones(size(TsMdlData,1),1), TsMdlData.Rwd,25,TsColors(AllTs==TSi,:).*ones(size(TsMdlData,1),1),'filled');
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            ylim([min(max(-1000, min(TsMdlData.Rwd)), min(ylim)) 0]);
            xlim([min(AllTs)-0.005, max(AllTs)+0.005]);
        end
        
    end
    
end

% Result Tables
TableMeansE = array2table(TableMeansE);
TableMeansE.Properties.VariableNames = {'ExpSetN', 'Ts', ...
    'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
    'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};
TableMediansE = array2table(TableMediansE);
TableMediansE.Properties.VariableNames = {'ExpSetN', 'Ts', ...
      'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
      'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};
TableWavgsE = array2table(TableWavgsE);
TableWavgsE.Properties.VariableNames = {'ExpSetN', 'Ts', ...
    'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
    'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};

%}


%% NTraj & Ts loop
%{
fig5 = figure(5);
fig5.Color = [1,1,1];
fig5.Units = 'normalized';
fig6 = figure(6);
fig6.Color = [1,1,1];
fig6.Units = 'normalized';
fig6.Position(3:4) = [0.4 0.35];

TableMeansT = [];
TableMediansT = [];
TableWavgsT = [];

spi = 1;

MdlData4 = MdlData(MdlData.nCOM==4,:);
for NTi=AllTrajs'
    
    %NTrajData = MdlData4(MdlData4.NExp==NTi,:);
    NTrajData = MdlData(MdlData.NExp==NTi,:);
    
    for TSi=AllTs'
        
        TsMdlData = NTrajData(NTrajData.Ts==TSi,:);
        ThTsMdlData = table2array(TsMdlData(:,ThCols))';
        RwTsMdlData = 1 - (TsMdlData.Rwd - min(TsMdlData.Rwd))/(range(TsMdlData.Rwd)+eps);
        PlotX = kron(ones(size(TsMdlData,1),1), (1:7)');
        PlotC = kron(RwTsMdlData, ones(7,1));
        
        if size(ThTsMdlData,2)>0
            dw = REPSupdate(TsMdlData.Rwd);
            
            ThTsMdlMean = mean(ThTsMdlData,2);
            ThTsMdlMedian = median(ThTsMdlData,2);
            ThTsMdlWavg = sum(ThTsMdlData.*dw,2)/sum(dw);
            
            TableMeansT = [TableMeansT; NTi TSi ThTsMdlMean'];
            TableMediansT = [TableMediansT; NTi TSi ThTsMdlMedian'];
            TableWavgsT = [TableWavgsT; NTi TSi ThTsMdlWavg'];
            
            figure(6);
            subplot(1,9,1:3)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            %scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            scatter((1:7)', ThTsMdlWavg, 50, 'om', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([0.5,3.5]);
            xticklabels({'$k_x$','$k_y$','$k_z$'});
            ylabel('Value [N/m]','Interpreter','latex','FontSize',10);
            subplot(1,9,5:7)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            %scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            scatter((1:7)', ThTsMdlWavg, 50, 'om', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([3.5,6.5]);
            xticklabels({'$b_x$','$b_y$','$b_z$'});
            ylabel('Value [Ns/m]','Interpreter','latex','FontSize',10);
            subplot(1,9,9)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            %scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            scatter((1:7)', ThTsMdlWavg, 50, 'om', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([6.5,7.5]);
            xticklabels({'','$\Delta l_{0z}$',''});
            ylabel('Value [m]','Interpreter','latex','FontSize',10);
            
            %{
            figure(5);
            subplot(length(AllCOMSizes), length(AllTs), spi)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            if size(ThTsMdlData,2)>1
                hold on
                boxplot(ThTsMdlData');
                hold off
            end
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            %}
            %{
            figure(5);
            hold on
            scatter3(NTi*ones(size(TsMdlData,1),1),TSi*ones(size(TsMdlData,1),1),TsMdlData.Rwd,25,RwTsMdlData,'filled');
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            zlim([max(-1000, min(zlim)) 0]);
            %}
            figure(5);
            subplot(1,2,1)
            hold on
            scatter(NTi*ones(size(TsMdlData,1),1), TsMdlData.Rwd,25,RwTsMdlData,'filled');
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            ylim([min(max(-1000, min(TsMdlData.Rwd)), min(ylim)) 0]);
            xlim([min(AllTrajs)-1, max(AllTrajs)+1]);
            subplot(1,2,2)
            hold on
            scatter(TSi*ones(size(TsMdlData,1),1), TsMdlData.Rwd,25,RwTsMdlData,'filled');
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            ylim([min(max(-1000, min(TsMdlData.Rwd)), min(ylim)) 0]);
            xlim([min(AllTs)-0.005, max(AllTs)+0.005]);
        end
        spi = spi+1;

    end
    
end

% Result Tables
TableMeansT = array2table(TableMeansT);
TableMeansT.Properties.VariableNames = {'NTraj', 'Ts', ...
  'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
  'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};
TableMediansT = array2table(TableMediansT);
TableMediansT.Properties.VariableNames = {'NTraj', 'Ts', ...
    'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
    'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};
TableWavgsT = array2table(TableWavgsT);
TableWavgsT.Properties.VariableNames = {'NTraj', 'Ts', ...
    'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
    'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};
%}


%% nCOM & Ts loop

fig2 = figure(2);
fig2.Color = [1,1,1];
fig2.Units = 'normalized';
fig2.Position = [0.05 0.1 0.9 0.35];
fig3 = figure(3);
fig3.Color = [1,1,1];
fig3.Units = 'normalized';
fig3.Position = [0.05 0.5 0.35 0.35];
fig4 = figure(4);
fig4.Color = [1,1,1];
fig4.Units = 'normalized';
fig4.Position = [0.55 0.5 0.4 0.35];

TableMeansN = [];
TableMediansN = [];
TableWavgsN = [];

cmj = jet(1000);
spi = 1;

for NCi=AllCOMSizes'
    
    COMSizeData = MdlData(MdlData.nCOM==NCi,:);
    
    for TSi=AllTs'
        
        TsMdlData = COMSizeData(COMSizeData.Ts==TSi,:);
        ThTsMdlData = table2array(TsMdlData(:,ThCols))';
        RwTsMdlData = 1 - (TsMdlData.Rwd - min(TsMdlData.Rwd))/(range(TsMdlData.Rwd)+eps);
        PlotX = kron(ones(size(TsMdlData,1),1), (1:7)');
        PlotC = kron(RwTsMdlData, ones(7,1));
        
        if size(ThTsMdlData,2)>0
            dw = REPSupdate(TsMdlData.Rwd);
            
            for NTi=AllTrajs'

                TrajRows = TsMdlData.NExp==NTi;
                TrajMdlData = TsMdlData(TrajRows,:);
                TrajTh = ThTsMdlData(:,TrajRows);
                TrajRw = RwTsMdlData(TrajRows);

                figure(2);
                
                ThYlabels={'$k_x$ Value [N/m]','$k_y$ Value [N/m]','$k_z$ Value [N/m]', ...
                        '$b_x$ Value [Ns/m]','$b_y$ Value [Ns/m]','$b_z$ Value [Ns/m]', ...
                        '$\Delta l_{0z}$ Value [m]'};
                for thtraj=1:7
                    subplot(1,7,thtraj)
                    if(NTi ~= AllTrajs(1)), hold on; end
                    scatter(NTi*ones(size(TrajMdlData,1),1), TrajTh(thtraj,:), 25, TrajRw, 'filled');
                    hold off
                    grid on;
                    box on;
                    set(gca, 'TickLabelInterpreter','latex');
                    xlim([min(AllTrajs)-0.5,max(AllTrajs)+0.5]);
                    %xticklabels({'$k_x$','$k_y$','$k_z$'});
                    ylabel(ThYlabels{thtraj},'Interpreter','latex','FontSize',10);
                end

            end
            colormap(cmj(100:900,:));
            cb2ticks = min(abs(TsMdlData.Rwd))+(0:0.1:1)*range(TsMdlData.Rwd);
            cb2ticks = round(cb2ticks*100)/100;
            cb2 = colorbar('Position',[0.925 0.125 0.02 0.8], ...
                  'TickLabels',cb2ticks, 'TickLabelInterpreter','latex');
            cb2.Label.String = 'Cost $=|$Reward$|$';
            cb2.Label.Interpreter = 'latex';

            
            ThTsMdlMean = mean(ThTsMdlData,2);
            ThTsMdlMedian = median(ThTsMdlData,2);
            ThTsMdlWavg = sum(ThTsMdlData.*dw,2)/sum(dw);
            
            TableMeansN = [TableMeansN; NCi TSi ThTsMdlMean'];
            TableMediansN = [TableMediansN; NCi TSi ThTsMdlMedian'];
            TableWavgsN = [TableWavgsN; NCi TSi ThTsMdlWavg'];
            
            figure(4);
            subplot(1,9,1:3)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            %scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            %scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            scatter((1:7)', ThTsMdlWavg, 50, 'om', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([0.5,3.5]);
            xticklabels({'$k_x$','$k_y$','$k_z$'});
            ylabel('Value [N/m]','Interpreter','latex','FontSize',10);
            subplot(1,9,5:7)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            %scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            %scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            scatter((1:7)', ThTsMdlWavg, 50, 'om', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([3.5,6.5]);
            xticklabels({'$b_x$','$b_y$','$b_z$'});
            ylabel('Value [Ns/m]','Interpreter','latex','FontSize',10);
            subplot(1,9,9)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            %scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            %scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            scatter((1:7)', ThTsMdlWavg, 50, 'om', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([6.5,7.5]);
            xticklabels({'','$\Delta l_{0z}$',''});
            ylabel('Value [m]','Interpreter','latex','FontSize',10);
            colormap(cmj(100:900,:));
            
            cb4ticks = min(abs(TsMdlData.Rwd))+(0:0.1:1)*range(TsMdlData.Rwd);
            cb4ticks = round(cb4ticks*100)/100;
            cb4 = colorbar('Position',[0.87 0.11 0.035 0.815], ...
                  'TickLabels',cb2ticks, 'TickLabelInterpreter','latex');
            cb4.Label.String = 'Cost $=|$Reward$|$';
            cb4.Label.Interpreter = 'latex';
            cb4.Label.FontSize=10;
            
            %{
            figure(3);
            subplot(length(AllCOMSizes), length(AllTs), spi)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            if size(ThTsMdlData,2)>1
                hold on
                boxplot(ThTsMdlData');
                hold off
            end
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            %}
            %{
            figure(3);
            hold on
            scatter3(NCi*ones(size(TsMdlData,1),1),TSi*ones(size(TsMdlData,1),1),TsMdlData.Rwd,25,RwTsMdlData,'filled');
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            zlim([-1000 0]);
            %}
            
            figure(3);
            subplot(1,2,1)
            if(spi ~= 1), hold on; end
            scatter(NCi*ones(size(TsMdlData,1),1), TsMdlData.Rwd,25,RwTsMdlData,'filled');
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            ylim([min(max(-1200, min(TsMdlData.Rwd)), min(ylim)) 0]);
            xlim([min(AllCOMSizes)-0.5, max(AllCOMSizes)+0.5]);
            xlabel('Model Side Size ($n$)', 'Interpreter','latex')
            ylabel('Reward', 'Interpreter','latex')
            subplot(1,2,2)
            if(spi ~= 1), hold on; end
            scatter(TSi*ones(size(TsMdlData,1),1), TsMdlData.Rwd,25,RwTsMdlData,'filled');
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            ylim([min(max(-1200, min(TsMdlData.Rwd)), min(ylim)) 0]);
            xlim([min(AllTs)-0.0025, max(AllTs)+0.0025]);
            xlabel('$T_s$ [s]', 'Interpreter','latex')
            ylabel('Reward', 'Interpreter','latex')
            %xticklabels({'10','15','20','25'});
            colormap(cmj(100:900,:));
        end
        spi = spi+1;

    end
    
end

% Result Tables
TableMeansN = array2table(TableMeansN);
TableMeansN.Properties.VariableNames = {'MdlSz', 'Ts', ...
  'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
  'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};
TableMediansN = array2table(TableMediansN);
TableMediansN.Properties.VariableNames = {'MdlSz', 'Ts', ...
    'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
    'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};
TableWavgsN = array2table(TableWavgsN);
TableWavgsN.Properties.VariableNames = {'MdlSz', 'Ts', ...
    'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
    'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};


%% Save final table
if(saveTable==1)
    writetable(TableMediansN,'ThetaMdl_LUT.csv');
    fprintf('Saved ThetaMdl_LUT.csv\n');
end

  
%% Linear Regression

AllResults = table2array(TableMediansN);
RX = AllResults(:,1:2);
RX = [RX, RX.^2 RX(:,1).*RX(:,2)];
RY = AllResults(:,3:end);

RK = (RX'*RX)\RX'*RY;

%{x
rSz = 7;
rTs = 0.025;
rX  = [rSz rTs rSz^2 rTs^2 rSz*rTs];
rY  = rX*RK;
%}








