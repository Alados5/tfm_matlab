clear; clc;

%% Initialization

MdlData = readtable('ThetaModelLUT.csv');
ThCols = contains(MdlData.Properties.VariableNames, 'Th_');
ThColNames = MdlData.Properties.VariableNames(ThCols);
MainColNames = {'ExpSetN','NExp','NTrial','Ts','nSOM','nCOM'};

% Filter unnecessary columns
MdlData = MdlData(:, [MainColNames(:)', ThColNames(:)']);
ThCols = contains(MdlData.Properties.VariableNames, 'Th_');

% Remove trajectories not fit for learning
MdlData = MdlData(MdlData.ExpSetN~=11,:);
MdlData = MdlData(~sum(MdlData.NExp==[1:3,6,10:13], 2), :);
MdlData = MdlData(~sum(MdlData.NTrial==3, 2), :);

AllExpSets = unique(MdlData.ExpSetN);
AllCOMSizes = unique(MdlData.nCOM);
AllTs = unique(MdlData.Ts);
TotalN1 = length(AllTs)*length(AllExpSets);
TotalN2 = length(AllTs)*length(AllCOMSizes);

TableMeansE = [];
TableMediansE = [];
TableMeansN = [];
TableMediansN = [];

%% ExpSetN & Ts loop

fig1 = figure(1);
fig1.Color = [1,1,1];
fig2 = figure(2);
fig2.Color = [1,1,1];
spi = 1;
colormap jet;

for ESi=AllExpSets'
    
    ExpSetData = MdlData(MdlData.ExpSetN==ESi,:);
    
    for TSi=AllTs'
        
        TsMdlData = ExpSetData(ExpSetData.Ts==TSi,:);
        ThTsMdlData = table2array(TsMdlData(:,ThCols))';
        PlotX = kron(ones(size(TsMdlData,1),1), (1:7)');
        PlotC = kron((1:size(TsMdlData,1))', ones(7,1));
        
        if size(ThTsMdlData,2)>0
            ThTsMdlMean = mean(ThTsMdlData,2);
            ThTsMdlMedian = median(ThTsMdlData,2);
            
            TableMeansE = [TableMeansE; ESi TSi ThTsMdlMean'];
            TableMediansE = [TableMediansE; ESi TSi ThTsMdlMedian'];
                        
            figure(2);
            subplot(1,2,1)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([1,3]);
            subplot(1,2,2)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([4,7]);
            
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
        end
        spi = spi+length(AllExpSets);

    end
    spi=mod(spi+1, TotalN1);
    
end


%% nCOM & Ts loop

fig3 = figure(3);
fig3.Color = [1,1,1];
fig4 = figure(4);
fig4.Color = [1,1,1];
spi = 1;

for NCi=AllCOMSizes'
    
    COMSizeData = MdlData(MdlData.nCOM==NCi,:);
    
    for TSi=AllTs'
        
        TsMdlData = COMSizeData(COMSizeData.Ts==TSi,:);
        ThTsMdlData = table2array(TsMdlData(:,ThCols))';
        PlotX = kron(ones(size(TsMdlData,1),1), (1:7)');
        PlotC = kron((1:size(TsMdlData,1))', ones(7,1));
        
        if size(ThTsMdlData,2)>0
            ThTsMdlMean = mean(ThTsMdlData,2);
            ThTsMdlMedian = median(ThTsMdlData,2);
            
            TableMeansN = [TableMeansN; NCi TSi ThTsMdlMean'];
            TableMediansN = [TableMediansN; NCi TSi ThTsMdlMedian'];
                        
            figure(4);
            subplot(1,2,1)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([1,3]);
            subplot(1,2,2)
            scatter(PlotX, ThTsMdlData(:), 25, PlotC, 'filled');
            hold on
            scatter((1:7)', ThTsMdlMean, 50, 'om', 'LineWidth',1);
            scatter((1:7)', ThTsMdlMedian, 50, 'xm', 'LineWidth',1);
            hold off
            grid on;
            box on;
            set(gca, 'TickLabelInterpreter','latex');
            xlim([4,7]);
            
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
        end
        spi = spi+1;

    end
    
end


%% Result Tables

TableMeansE = array2table(TableMeansE);
TableMeansE.Properties.VariableNames = {'ExpSetN', 'Ts', ...
    'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
    'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};
TableMediansE = array2table(TableMediansE);
TableMediansE.Properties.VariableNames = {'ExpSetN', 'Ts', ...
      'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
      'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};
  
  
TableMeansN = array2table(TableMeansN);
TableMeansN.Properties.VariableNames = {'MdlSz', 'Ts', ...
  'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
  'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};
TableMediansN = array2table(TableMediansN);
TableMediansN.Properties.VariableNames = {'MdlSz', 'Ts', ...
    'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
    'Th_damping_x','Th_damping_y','Th_damping_z', 'Th_z_sum'};

% Save final table
writetable(TableMediansN,'LearntModelParams.csv');
fprintf('Saved LearntModelParams.csv\n');

  
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








