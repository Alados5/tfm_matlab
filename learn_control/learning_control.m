close all; clc; clear;

%% Initialization
ExpSetN = 1;
SimType = 'LIN'; %LIN, NL, RTM
ExpSetNote = '';
NTraj = 6;
Ts = 0.025;
Hp = 20;
Wv = 0.3;
sigmaD = 0.0;
sigmaN = 0.0;
nSOM = 4;
nCOM = 4;
nNLM = 10;

SOM_ThetaExp = [4,8,2];
COM_ThetaExp = [4,8,2];

e0 = 60;
minRwd = -100;
NSamples = 10;
NEpochs = 5;
UseLambda = 1;

Plot3DTraj = 0;
Plot2DTraj = 1;

% Parameters to learn and defaults:
%  [ubound=5*1e-3, gbound=0*1e-3, W_Q=1, W_T=1, W_R=10]
ThMask = [0.001 0.001 1 1 1]';
ThW = 3:5;


% Load parameter table and select corresponding row
ThetaModelLUT = readtable('./ThetaModelLUT.csv');
LUT_SOM_id = (ThetaModelLUT.ExpSetN == SOM_ThetaExp(1)) & ...
             (ThetaModelLUT.NExp == SOM_ThetaExp(2)) & ...
             (ThetaModelLUT.NTrial == SOM_ThetaExp(3)) & ... 
             (ThetaModelLUT.Ts == Ts) & (ThetaModelLUT.nCOM == nSOM);
LUT_COM_id = (ThetaModelLUT.ExpSetN == COM_ThetaExp(1)) & ...
             (ThetaModelLUT.NExp == COM_ThetaExp(2)) & ...
             (ThetaModelLUT.NTrial == COM_ThetaExp(3)) & ...
             (ThetaModelLUT.Ts == Ts) & (ThetaModelLUT.nCOM == nCOM);
LUT_SOM = ThetaModelLUT(LUT_COM_id, :);
LUT_COM = ThetaModelLUT(LUT_COM_id, :);
if (size(LUT_COM,1) > 1 || size(LUT_SOM,1) > 1)
    error("There are multiple rows with same experiment parameters.");
elseif (size(LUT_COM,1) < 1 || size(LUT_SOM,1) < 1)
    error("There are no saved experiments with those parameters.");
else
    paramsSOM = table2array(LUT_SOM(:, contains(LUT_SOM.Properties.VariableNames, 'Th_')));
    paramsCOM = table2array(LUT_COM(:, contains(LUT_COM.Properties.VariableNames, 'Th_')));
end


% Options
opts = struct();
opts.NTraj = NTraj;
opts.Ts = Ts;
opts.Hp = Hp;
opts.Wv = Wv;
opts.sigmaD = sigmaD;
opts.sigmaN = sigmaN;
opts.nSOM = nSOM;
opts.nCOM = nCOM;
opts.nNLM = nNLM;
opts.paramsSOM = paramsSOM;
opts.paramsCOM = paramsCOM;
opts.xbound = 1.5;
% ----------


if strcmp(SimType, 'LIN')
    SimTypeN = 0;
    wvname = '';
elseif strcmp(SimType, 'NL')
    SimTypeN = 1;
    wvname = '';
elseif strcmp(SimType, 'RTM') || strcmp(SimType, 'RT')
    SimTypeN = 2;
    wvname = ['_wv', num2str(Wv*100)];
else
    error(['Simulation type "',SimType,'" is not a valid option.']);
end

dirname = ['Exps',num2str(ExpSetN), '_',SimType, ExpSetNote, ...
           '/traj',num2str(NTraj),'_ts',num2str(Ts*1000), ...
           '_hp',num2str(Hp), wvname, '_ns',num2str(nSOM), '_nc',num2str(nCOM)];


if e0==0
    % Initial seed [xbound; ubound; gbound; WQ; WT; WR]
    mw0 = [5; 0;          0.5; 0.5; 0.5];
    Sw0 = diag([2; 0.1;   0.5; 0.5; 0.5]);
else
    prevrange = [num2str(e0-NEpochs),'-',num2str(e0)];

    MWans = load([dirname,'/MW_',prevrange,'.mat']);
    MWans = MWans.MW;
    mw0 = MWans(:,:,end);
    SWans = load([dirname,'/SW_',prevrange,'.mat']);
    SWans = SWans.SW;
    Sw0 = SWans(:,:,end);
    RWans = load([dirname,'/RW_',prevrange,'.mat']);
    RWans = RWans.RW;
    RWprev = RWans(:,:,end);
    THans = load([dirname,'/TH_',prevrange,'.mat']);
    THans = THans.TH;
    THprev = THans(:,:,end);
end

NParams = length(mw0);

MW = zeros(NParams,1,NEpochs+1);
SW = zeros(NParams,NParams,NEpochs+1);
TH = zeros(NParams,NSamples,NEpochs);
RW = zeros(1,NSamples,NEpochs);
XR = cell(1,NSamples,NEpochs);

MW(:,:,1) = mw0;
SW(:,:,1) = Sw0;


%% Main Learning Loop
mw = mw0;
Sw = Sw0;
epoch = 1;
while epoch <= NEpochs
    
    wghts_ep = zeros(NParams,NSamples);
    rwrds_ep = zeros(1,NSamples);

    fprintf(['\nEpoch: ' num2str(e0+epoch),'\n-----------------------']);
    
    for i=1:NSamples
        
        fprintf(['\nEpoch ', num2str(e0+epoch), ' - Sample: ' num2str(i),'\n']);
        
        lambda = mean(svd(Sw))/5;
        SwL = Sw + UseLambda*eye(NParams)*lambda;
        %SwL(7,7) = Sw(7,7)*1.1;
        
        thetai = mvnrnd(mw, SwL)';
        th_maxW = max(thetai(ThW));
        thetai(ThW) = thetai(ThW)/th_maxW;
        while any(thetai<0) || any(thetai(ThW)>1) || all(thetai(ThW)==0)
            thetai = mvnrnd(mw, SwL)';
            th_maxW = max(thetai(ThW));
            thetai(ThW) = thetai(ThW)/th_maxW;
        end

        theta = thetai.*ThMask;
        fprintf([' Theta: [',num2str(theta',5),']\n']);
        
        if SimTypeN==2
            [Rwd, AllData] = simulation_cl_rtm_theta(theta, opts);
        elseif SimTypeN==1
            [Rwd, AllData] = simulation_cl_nl_theta(theta, opts);
        else
            [Rwd, AllData] = simulation_cl_lin_theta(theta, opts);
        end
        Rwd = max(Rwd, minRwd);

        wghts_ep(:,i) = thetai;
        rwrds_ep(i) = Rwd;
        XR{1,i,epoch} = AllData;
    end
    
    if isequal(unique(rwrds_ep), minRwd)
        fprintf('\nEPOCH HAD NO SUCCESSFUL RESULTS. RE-DOING SAME EPOCH \n');
        pause(1);
    else
        if epoch==1
            if e0==0
                dw = REPSupdate([rwrds_ep, rwrds_ep]);
                weights2 = [wghts_ep wghts_ep];
            else
                dw = REPSupdate([rwrds_ep, RWprev]);
                weights2 = [wghts_ep THprev];
            end
        else
            dw = REPSupdate([rwrds_ep, RW(:,:,epoch-1)]);
            weights2 = [wghts_ep TH(:,:,epoch-1)];
        end
        
        Z = (sum(dw)*sum(dw) - sum(dw .^ 2))/sum(dw);
        mw = sum(bsxfun(@times, weights2', dw),1)'./sum(dw);
        mw(ThW) = mw(ThW)/max(mw(ThW));
        summ = 0;
        for ak = 1:size(weights2,2)
            summ = summ + dw(ak)*((weights2(:,ak)-mw)*(weights2(:,ak)-mw)');%/sum(dw);
        end
        Sw = summ./(Z+1e-9);

        MW(:,:,epoch+1) = mw;
        SW(:,:,epoch+1) = Sw;
        TH(:,:,epoch)   = wghts_ep;
        RW(:,:,epoch)   = rwrds_ep;
        
        epoch = epoch+1;
    end

end

% Save data
if ~isfolder(dirname)
    mkdir(dirname);
end

epochrange = [num2str(e0),'-',num2str(e0+NEpochs)];

save([dirname,'/MW_',epochrange,'.mat'],'MW');
save([dirname,'/SW_',epochrange,'.mat'],'SW');
save([dirname,'/RW_',epochrange,'.mat'],'RW');
save([dirname,'/TH_',epochrange,'.mat'],'TH');


%% Execution of mean weights of each epoch
MW2D = permute(MW,[1,3,2]);
RWMW = zeros(1, size(MW2D,2));

fprintf('\nExecuting resulting means per epoch...\n-----------------------');

for epoch=1:size(MW2D,2)
    fprintf(['\nEpoch: ', num2str(e0+epoch-1), '\t|']);
    theta = (MW2D(:,epoch).*ThMask)';
    if SimTypeN==2
        [Rwd, AllData] = simulation_cl_rtm_theta(theta, opts);
    elseif SimTypeN==1
        [Rwd, AllData] = simulation_cl_nl_theta(theta, opts);
    else
        [Rwd, AllData] = simulation_cl_lin_theta(theta, opts);
    end
    RWMW(epoch) = Rwd;
end
ThLearnt = (MW2D(:,end).*ThMask)';
fprintf(['\nLearnt Theta: [',num2str(ThLearnt,5),']\n']);

% Save data
save([dirname,'/RWMW_',epochrange,'.mat'],'RWMW');
fprintf('Saved data files. \n');


%% Reward Plots

hf1 = figure(1);
hf1.Color = [1,1,1];
hf1.Units = 'normalized';
hf1.Position = [0 0.5 0.5 0.4];

% Evolution of sample mean rewards
RW2D = permute(RW, [2,3,1]);
RWm = mean(RW2D);

subplot(2,1,1);
plot(e0:e0+NEpochs-1, RWm, 'LineWidth',1);
grid on
set(gca,'TickLabelInterpreter','latex');
xlim([e0 e0+NEpochs]);
xlabel('Epoch', 'Interpreter','latex');
ylabel('Reward', 'Interpreter','latex');
title('\textbf{Mean of sample rewards}', 'Interpreter','latex');

% Evolution of resulting means

subplot(2,1,2);
plot(e0:e0+NEpochs, RWMW, 'LineWidth',1);
hold on
plot([e0 e0+NEpochs], [RWMW(1) RWMW(1)], '--k');
hold off
grid on
set(gca,'TickLabelInterpreter','latex');
xlabel('Epoch', 'Interpreter','latex');
ylabel('Reward', 'Interpreter','latex');
title('\textbf{Reward of resulting mean}', 'Interpreter','latex');
if (min(min(hf1.Children.YLim)) < minRwd)
    ylim([minRwd, min(max(RWMW)+10,0)]);
end


% Full experiment
hf2 = figure(2);
hf2.Color = [1,1,1];
hf2.Units = 'normalized';
hf2.Position = [0 0.05 0.5 0.4];

full_exp_plot(dirname, 2);
if (min(min(hf2.Children.YLim)) < minRwd)
    ylim([minRwd, min(max(RWMW)+10,0)]);
end


%% 3D & 2D Plots
tSim = 0:Ts:size(AllData.xSOM,2)*Ts-Ts;
SOM_node_ctrl = [opts.nSOM*(opts.nSOM-1)+1, opts.nSOM*opts.nSOM];
COM_node_ctrl = [opts.nCOM*(opts.nCOM-1)+1, opts.nCOM*opts.nCOM];
NLM_node_ctrl = [opts.nNLM*(opts.nNLM-1)+1, opts.nNLM*opts.nNLM];
coord_ctrlS = [SOM_node_ctrl SOM_node_ctrl+opts.nSOM^2 SOM_node_ctrl+2*opts.nSOM^2];
coord_ctrlC = [COM_node_ctrl COM_node_ctrl+opts.nCOM^2 COM_node_ctrl+2*opts.nCOM^2];
coord_ctrlN = [NLM_node_ctrl NLM_node_ctrl+opts.nNLM^2 NLM_node_ctrl+2*opts.nNLM^2];
coord_lcS = [1 opts.nSOM 1+opts.nSOM^2 opts.nSOM^2+opts.nSOM 2*opts.nSOM^2+1 2*opts.nSOM^2+opts.nSOM];
coord_lcC = [1 opts.nCOM 1+opts.nCOM^2 opts.nCOM^2+opts.nCOM 2*opts.nCOM^2+1 2*opts.nCOM^2+opts.nCOM]; 
coord_lcN = [1 opts.nNLM 1+opts.nNLM^2 opts.nNLM^2+opts.nNLM 2*opts.nNLM^2+1 2*opts.nNLM^2+opts.nNLM]; 


Ref_l = load(['../data/trajectories/ref_',num2str(NTraj),'L.csv']);
Ref_r = load(['../data/trajectories/ref_',num2str(NTraj),'R.csv']);
nPtRef = size(Ref_l,1);


% 3D trajectories (uc & lc)
if(Plot3DTraj==1)
    hf3 = figure(3);
    hf3.Color = [1,1,1];

    plot3(AllData.xSOM(coord_ctrlS(1:2),:)', AllData.xSOM(coord_ctrlS(3:4),:)', AllData.xSOM(coord_ctrlS(5:6),:)');
    hold on
    plot3(AllData.xSOM(coord_lcS(1:2),:)', AllData.xSOM(coord_lcS(3:4),:)', AllData.xSOM(coord_lcS(5:6),:)');
    plot3(Ref_l(:,1), Ref_l(:,2), Ref_l(:,3), '--k', 'linewidth',1);
    plot3(Ref_r(:,1), Ref_r(:,2), Ref_r(:,3), '--k', 'linewidth',1);
    hold off
    box on; grid on;
    axis equal
    set(gca,'TickLabelInterpreter','latex');
end

% 2D evolutions
if(Plot2DTraj==1)
    hf4 = figure(4);
    hf4.Color = [1,1,1];
    hf4.Units = 'normalized';
    hf4.Position = [0.5 0.05 0.5 0.65];
    
    subplot(15,2,1:2:12);
    plot(tSim, AllData.xSOM(coord_ctrlS([1,3,5]),:), 'linewidth',1.5);
    grid on
    set(gca,'TickLabelInterpreter','latex');
    ylim41 = ylim;
    
    subplot(15,2,2:2:12);
    plot(tSim, AllData.xSOM(coord_ctrlS([2,4,6]),:), 'linewidth',1.5);
    grid on
    set(gca,'TickLabelInterpreter','latex');
    ylim42 = ylim;
    
    ylim4high = [min(ylim41(1), ylim42(1)), max(ylim41(2), ylim42(2))];
    subplot(15,2,1:2:12);
    ylim(ylim4high);
    subplot(15,2,2:2:12);
    ylim(ylim4high);

    subplot(15,2,17:2:28);
    plot(tSim, AllData.xSOM(coord_lcS([1,3,5]),:), 'linewidth',1.5);
    hold on
    plot(tSim, Ref_l, '--k', 'linewidth',1);
    hold off
    grid on
    set(gca,'TickLabelInterpreter','latex');
    ylim43 = ylim;

    subplot(15,2,18:2:28);
    pa4som = plot(tSim, AllData.xSOM(coord_lcS([2,4,6]),:), 'linewidth',1.5);
    hold on
    pa4ref = plot(tSim, Ref_r, '--k', 'linewidth',1);
    hold off
    grid on
    set(gca,'TickLabelInterpreter','latex');
    ylim44 = ylim;

    ylim4low = [min(ylim43(1), ylim44(1)), max(ylim43(2), ylim44(2))];
    subplot(15,2,17:2:28);
    ylim(ylim4low);
    subplot(15,2,18:2:28);
    ylim(ylim4low);

    Lgnd4 = legend([pa4som; pa4ref(1)], ...
                   '$x_{SOM}$','$y_{SOM}$','$z_{SOM}$', 'Ref', ...
                   'Orientation', 'horizontal', 'Interpreter', 'latex');
    Lgnd4.Position(1) = 0.5-Lgnd4.Position(3)/2;
    Lgnd4.Position(2) = 0.05;

end


%% Update LUT for Parameters
ThetaCtrlLUT = readtable('ThetaControlLUT.csv');
LUT_Exp_id = (ThetaCtrlLUT.ExpSetN == ExpSetN) & ...
             (ThetaCtrlLUT.NTraj == NTraj) & ...
             (ThetaCtrlLUT.SimType == SimTypeN) & ...
             (ThetaCtrlLUT.Ts == Ts) & (ThetaCtrlLUT.Hp == Hp) & ...
             (ThetaCtrlLUT.nSOM == nSOM) & (ThetaCtrlLUT.nCOM == nCOM) & ...
             (ThetaCtrlLUT.NEpochs == NEpochs) & (ThetaCtrlLUT.NSamples == NSamples) & ...
             (ThetaCtrlLUT.MinRwd == minRwd);
LUT_Exp = ThetaCtrlLUT(LUT_Exp_id, :);

if (size(LUT_Exp,1) > 1)
    error("There are multiple rows with same experiment parameters.");
elseif (size(LUT_Exp,1) == 1)
    % Update experiment row
    LUT_row = find(LUT_Exp_id);
    ThetaCtrlLUT(LUT_row,'LastEpoch') = {e0+NEpochs};
    ThetaCtrlLUT(LUT_row,'Th_ubound') = {ThLearnt(1)};
    ThetaCtrlLUT(LUT_row,'Th_gbound') = {ThLearnt(2)};
    ThetaCtrlLUT(LUT_row,'Th_WQ') = {ThLearnt(3)};
    ThetaCtrlLUT(LUT_row,'Th_WT') = {ThLearnt(4)};
    ThetaCtrlLUT(LUT_row,'Th_WR') = {ThLearnt(5)};
else
    % Add experiment row
    LUT_row = size(ThetaCtrlLUT,1)+1;
    ThetaCtrlLUT(LUT_row,'ExpSetN') = {ExpSetN};
    ThetaCtrlLUT(LUT_row,'SimType') = {SimTypeN};
    ThetaCtrlLUT(LUT_row,'NTraj') = {NTraj};
    ThetaCtrlLUT(LUT_row,'Ts') = {Ts};
    ThetaCtrlLUT(LUT_row,'Hp') = {Hp};
    ThetaCtrlLUT(LUT_row,'nSOM') = {nSOM};
    ThetaCtrlLUT(LUT_row,'nCOM') = {nCOM};
    ThetaCtrlLUT(LUT_row,'LastEpoch') = {e0+NEpochs};
    ThetaCtrlLUT(LUT_row,'NEpochs') = {NEpochs};
    ThetaCtrlLUT(LUT_row,'NSamples') = {NSamples};
    ThetaCtrlLUT(LUT_row,'MinRwd') = {minRwd};
    
    ThetaCtrlLUT(LUT_row,'Th_ubound') = {ThLearnt(1)};
    ThetaCtrlLUT(LUT_row,'Th_gbound') = {ThLearnt(2)};
    ThetaCtrlLUT(LUT_row,'Th_WQ') = {ThLearnt(3)};
    ThetaCtrlLUT(LUT_row,'Th_WT') = {ThLearnt(4)};
    ThetaCtrlLUT(LUT_row,'Th_WR') = {ThLearnt(5)};
end
writetable(ThetaCtrlLUT,'ThetaControlLUT.csv');
fprintf('Updated ThetaControlLUT.csv\n');

