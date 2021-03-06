%{
Relative Entropy Policy Search algorithm applied to learn the optimal
weights of a Model Predictive Controller
- Original linear cloth model by David Parent, modified by Adri? Luque
- Original REPSupdate function by Adri? Colom?
- MPC and simulations by Adri? Luque
- Main code and other functions by Adri? Luque
- Last Updated: September 2021
%}

close all; clc; clear;


%% Initialization

% Learning options
ExpSet = 4;
SimType = 'LIN'; %LIN, NL, RTM
ExpNote = '_Det';
NTraj = 16; % Exps4: 6, 3, 13, 12, 16
Ts = 0.020;
Hp = 25;
Wv = 0.3;
nSOM = 4;
nCOM = 4;
nNLM = 10;
sigmaD = 0.0;
sigmaN = 0.0;
ubound  = 50*1e-3;  % (Enough Displ.)
gbound  = 0;        % (Eq. Constraint)

% MPC Structure options
opt_Du  = 1;  % 0=u,      1=Du
opt_Qa  = 0;  % 0=Qk,     1=Qa*Qk
opt_Rwd = 1;  % 1=RMSE,   2=Tov,           3=RMSE+Tov
opt_Wgh = 1;  % 1=[q r],  2=[qx qy qz r],  3=[qx qy qz k]

e0 = 0;
minRwd = -10;
NSamples = 100;
NEpochs = 2;
UseLambda = 1;

% Plotting options
Plot3DTraj = 0;
Plot2DTraj = (e0==0);


% Load parameter table and select corresponding row
ThetaModelLUT = readtable('../learn_model/ThetaMdl_LUT.csv');
LUT_SOM_id = (ThetaModelLUT.Ts == Ts) & (ThetaModelLUT.MdlSz == nSOM);
LUT_COM_id = (ThetaModelLUT.Ts == Ts) & (ThetaModelLUT.MdlSz == nCOM);
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
opts.nSOM = nSOM;
opts.nCOM = nCOM;
opts.nNLM = nNLM;
opts.sigmaD = sigmaD;
opts.sigmaN = sigmaN;
opts.ubound = ubound;
opts.gbound = gbound;
opts.opt_du = opt_Du;
opts.opt_Qa = opt_Qa;
opts.opt_Rwd = opt_Rwd;
opts.paramsSOM = paramsSOM;
opts.paramsCOM = paramsCOM;
% -------------------------


% Simulation type categorization
if strcmp(SimType, 'LIN')
    SimTypeN = 0;
    Wv=0; nNLM=0;
    wvname = '';
elseif strcmp(SimType, 'NL')
    SimTypeN = 1;
    Wv=0; nNLM=0;
    wvname = '';
elseif strcmp(SimType, 'RTM') || strcmp(SimType, 'RT')
    SimTypeN = 2;
    wvname = ['_wv', num2str(Wv*100)];
else
    error(['Simulation type "',SimType,'" is not a valid option.']);
end

% Weight structure categorization
if opt_Wgh==3
    %[qx; qy; qz; k]
    ThW = 1:3;
    ThM6 = @(mwi) [mwi(1) mwi(1) mwi(2) mwi(2) mwi(3) mwi(3);
                   mwi([1,1,2,2,3,3])'.*10^mwi(4)];
elseif opt_Wgh==2
    %[qx; qy; qz; r]
    ThW = 1:4;
    ThM6 = @(mwi) [mwi(1) mwi(1) mwi(2) mwi(2) mwi(3) mwi(3);
                   mwi(4) mwi(4) mwi(4) mwi(4) mwi(4) mwi(4)];
else
    %[q; r]
    ThW = 1:2;
    ThM6 = @(mwi) [mwi(1) mwi(1) mwi(1) mwi(1) mwi(1) mwi(1);
                   mwi(2) mwi(2) mwi(2) mwi(2) mwi(2) mwi(2)];
end

% Directory to save/load experiments
dirname = ['Exps',num2str(ExpSet), '/',num2str(SimTypeN), ...
           '_', SimType, ExpNote, ...
           '/traj',num2str(NTraj), '_ts',num2str(Ts*1000), ...
           '_hp',num2str(Hp), wvname, ...
           '_ns',num2str(nSOM), '_nc',num2str(nCOM)];


if e0==0
    % Initial seed
    if opt_Wgh==3
        %[qx; qy; qz; k]
        mw0 = [0.5; 0.5; 0.5; 0];
        Sw0 = diag([0.5; 0.5; 0.5; 1]);
    elseif opt_Wgh==2
        %[qx; qy; qz; r]
        mw0 = [0.5; 0.5; 0.5; 0.5];
        Sw0 = diag([0.5; 0.5; 0.5; 0.5]);
    else
        %[q; r]
        mw0 = [0.5; 0.5];
        Sw0 = diag([0.5; 0.5]);
    end
else
    % Load previous set of data (continue experiment)
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

% Initialize storage variables
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
        
        % Encourage exploration increasing variance
        lambda = mean(svd(Sw))/5;
        SwL = Sw + UseLambda*eye(NParams)*lambda;
        
        % Normalize weights so max=1, retry until valid sample
        thetai = mvnrnd(mw, SwL)';
        th_maxW = max(thetai(ThW));
        thetai(ThW) = thetai(ThW)/th_maxW;
        while any(thetai(ThW)<0) || any(thetai(ThW)>1) || all(thetai(ThW)==0)
            thetai = mvnrnd(mw, SwL)';
            th_maxW = max(thetai(ThW));
            thetai(ThW) = thetai(ThW)/th_maxW;
        end
        
        % Convert to 6 values for xxyyzz for both Q and R
        theta6 = ThM6(thetai);
        theta = struct;
        theta.Q = diag(theta6(1,:));
        theta.R = diag(theta6(2,:));
        fprintf([' Theta: [',num2str(thetai',5),']\n']);
        
        % Execute simulation
        if SimTypeN==2
            [Rwd, AllData] = sim_cl_rtm_theta(theta, opts);
        elseif SimTypeN==1
            [Rwd, AllData] = sim_cl_nl_theta(theta, opts);
        else
            [Rwd, AllData] = sim_cl_lin_theta(theta, opts);
        end
        
        % Saturate to minimum
        Rwd = max(Rwd, minRwd);

        % Save results
        wghts_ep(:,i) = thetai;
        rwrds_ep(i) = Rwd;
        XR{1,i,epoch} = AllData;
    end
    
    % Repeat epoch if no successful results (above saturation)
    if isequal(unique(rwrds_ep), minRwd)
        fprintf('\nEPOCH HAD NO SUCCESSFUL RESULTS. RE-DOING SAME EPOCH \n');
        pause(1);
        
    else
        % Add results of previous epoch to obtain relative weights
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
        
        % Update weights
        mw = sum(weights2'.*dw)'/sum(dw);
        mw(ThW) = mw(ThW)/max(mw(ThW));
        Z = (sum(dw)*sum(dw) - sum(dw .^ 2))/sum(dw);
        summ = 0;
        for ak = 1:size(weights2,2)
            summ = summ + dw(ak)*((weights2(:,ak)-mw)*(weights2(:,ak)-mw)');%/sum(dw);
        end
        Sw = summ./(Z+1e-9);

        % Save epoch results
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

fprintf(['\nExecuting resulting means per epoch...\n', ...
         '----------------------------------------']);

for epoch=1:size(MW2D,2)
    fprintf(['\nEpoch: ', num2str(e0+epoch-1), '\n']);
    
    theta6 = ThM6(MW2D(:,epoch));
    theta = struct;
    theta.Q = diag(theta6(1,:));
    theta.R = diag(theta6(2,:));

    if SimTypeN==2
        [Rwd, AllData] = sim_cl_rtm_theta(theta, opts);
    elseif SimTypeN==1
        [Rwd, AllData] = sim_cl_nl_theta(theta, opts);
    else
        [Rwd, AllData] = sim_cl_lin_theta(theta, opts);
    end
    RWMW(epoch) = Rwd;
end
ResLearnt = [AllData.eRMSE, AllData.eTov];
ThLearnt  = MW2D(:,end);
ThLearnt6 = ThM6(ThLearnt);
ThLearnt6 = [ThLearnt6(1,[1,3,5]) ThLearnt6(2,[1,3,5])];
fprintf(['\nLearnt Theta:  [',num2str(ThLearnt6,5),']']);
fprintf(['\nFinal Results: [RMSE = ',num2str(ResLearnt(1),5),', \t' ...
           'TOV = ',num2str(ResLearnt(2),5),']\n']);

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

% Distribution plots
%{
fooN = 100;
foo0 = abs(mvnrnd(mw0,Sw0,fooN));
foo1 = abs(mvnrnd(mw,Sw,fooN));
foo0n = foo0./max(foo0,[],2);
foo1n = foo1./max(foo1,[],2);

% Two weights
scatter(foo1n(:,1),foo1n(:,2),'ob','filled')
hold on
%scatter(foo0n(:,1),foo0n(:,2),'*r')
%scatter(foo0(:,1),foo0(:,2),'.m')
%scatter(foo1(:,1),foo1(:,2),'.c')
plot([1,1,0],[0,1,1],'--k')
hold off
grid on; box on; axis equal;
xlim([0 1.2])
ylim([0 1.2])

% XYZ Weights
scatter3(foo1n(:,1),foo1n(:,2),foo1n(:,3),'ob','filled')
hold on
%scatter3(foo0n(:,1),foo0n(:,2),foo0n(:,3),'*r')
plot3([1;1;0;0;1;1],[0;1;1;0;0;0],[1;1;1;1;1;0],'--k');
plot3([0;0;1;1;0;0],[1;0;0;1;1;1],[0;0;0;0;0;1],'--k');
plot3([0;0],[0;0],[1;0],'--k');
plot3([1;1],[1;1],[1;0],'--k');
hold off
grid on; box on; axis equal;
xlim([0 1.1])
ylim([0 1.1])
zlim([0 1.1])
xlabel('Qx');
ylabel('Qy');
zlabel('Qz');
%}


%% Update LUT for Parameters
DataTable = readtable('LearntCtrl_Data.csv');
Tbl_Exp_id = (DataTable.Exps == ExpSet) & ...
             (DataTable.Sim == SimTypeN) & ...
             (DataTable.NTraj == NTraj) & ...
             (DataTable.Ts == Ts) & (DataTable.Hp == Hp) & ...
             (DataTable.nSOM == nSOM) & (DataTable.nCOM == nCOM) & ...
             (DataTable.Wv == Wv) & (DataTable.nNLM == nNLM) & ...
             (DataTable.mRw == minRwd) & ...
             (DataTable.sD == sigmaD) & (DataTable.sN == sigmaN) & ...
             (DataTable.Du == opt_Du) & (DataTable.Qa == opt_Qa) & ...
             (DataTable.fRw == opt_Rwd) & (DataTable.Wgh == opt_Wgh);
Tbl_Exp = DataTable(Tbl_Exp_id, :);

if (size(Tbl_Exp,1) > 1)
    error("There are multiple rows with same experiment parameters.");
elseif (size(Tbl_Exp,1) == 1)
    % Update experiment row
    Tbl_row = find(Tbl_Exp_id);
    DataTable(Tbl_row,'nEps') = {e0+NEpochs};
    DataTable(Tbl_row,'Qx')   = {ThLearnt6(1)};
    DataTable(Tbl_row,'Qy')   = {ThLearnt6(2)};
    DataTable(Tbl_row,'Qz')   = {ThLearnt6(3)};
    DataTable(Tbl_row,'Rx')   = {ThLearnt6(4)};
    DataTable(Tbl_row,'Ry')   = {ThLearnt6(5)};
    DataTable(Tbl_row,'Rz')   = {ThLearnt6(6)};
    DataTable(Tbl_row,'RMSE') = {ResLearnt(1)};
    DataTable(Tbl_row,'TOV')  = {ResLearnt(2)};
else
    % Add experiment row
    Tbl_row = size(DataTable,1)+1;
    DataTable(Tbl_row,'Exps')  = {ExpSet};
    DataTable(Tbl_row,'Sim')   = {SimTypeN};
    DataTable(Tbl_row,'NTraj') = {NTraj};
    DataTable(Tbl_row,'Ts')    = {Ts};
    DataTable(Tbl_row,'Hp')    = {Hp};
    DataTable(Tbl_row,'Wv')    = {Wv};
    DataTable(Tbl_row,'nSOM')  = {nSOM};
    DataTable(Tbl_row,'nCOM')  = {nCOM};
    DataTable(Tbl_row,'nNLM')  = {nNLM};
    DataTable(Tbl_row,'sD')    = {sigmaD};
    DataTable(Tbl_row,'sN')    = {sigmaN};
    DataTable(Tbl_row,'Du')    = {opt_Du};
    DataTable(Tbl_row,'Qa')    = {opt_Qa};
    DataTable(Tbl_row,'fRw')   = {opt_Rwd};
    DataTable(Tbl_row,'Wgh')   = {opt_Wgh};
    DataTable(Tbl_row,'mRw')   = {minRwd};
    
    DataTable(Tbl_row,'nEps') = {e0+NEpochs};
    DataTable(Tbl_row,'Qx')   = {ThLearnt6(1)};
    DataTable(Tbl_row,'Qy')   = {ThLearnt6(2)};
    DataTable(Tbl_row,'Qz')   = {ThLearnt6(3)};
    DataTable(Tbl_row,'Rx')   = {ThLearnt6(4)};
    DataTable(Tbl_row,'Ry')   = {ThLearnt6(5)};
    DataTable(Tbl_row,'Rz')   = {ThLearnt6(6)};
    DataTable(Tbl_row,'RMSE') = {ResLearnt(1)};
    DataTable(Tbl_row,'TOV')  = {ResLearnt(2)};
end
writetable(DataTable,'LearntCtrl_Data.csv');
fprintf('Updated LearntCtrl_Data.csv\n');

