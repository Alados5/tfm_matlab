close all; clc; clear;

ExpSet = '11_Real_Cloth';
ExpDate = '2021_07_05';
NExp = 30;
NTrial = 1;
SimType = 'OL';
minRwd = -1000;
SizeSOM = 10;
SizeCOM = 4;
Ts = 0.015;

e0 = 20;
NSamples = 50;
NEpochs = 20;

Plot3DTraj = 0;
Plot2DTraj = 1;

ExpSetN = str2double(ExpSet(1:strfind(ExpSet,'_')-1));

% OL Options
opts = struct();
opts.nCOM = SizeCOM;
opts.Ts = Ts;
opts.ExpSetN = ExpSetN;
opts.NExp = NExp;
opts.NTrial = NTrial;
% NL SOM Model options
opts.nSOM = SizeSOM;
% Real Data options
opts.ExpDate = ExpDate;
% ----------


%{
% Original optimized for a 4x4
theta.stiffness = [-305.6028  -13.4221 -225.8987]; 
theta.damping = [-4.0042   -2.5735   -3.9090]; 
theta.z_sum = 0.0312;
%}


ThMask = [100 100 100 1 1 1 0.01]';
dirname = ['Exps',ExpSet,'/Exp',num2str(NExp),'_',num2str(NTrial),'_',num2str(Ts*1000),'ms'];
UseLambda = 1;
if e0==0
    % Initial seed for different sample times
    if (SizeCOM == 4 && Ts == 0.010)
        mw0 = [-2.5; -2.5; -2.5;    -4; -4; -4;   3];
        Sw0 = diag([0.5; 0.5; 0.5;   1; 1; 1;    1]);
        
    elseif (SizeCOM == 4 && Ts == 0.015)
        mw0 = [-2.0; -1.5; -2.0;     -3; -2.5; -3;     4];
        Sw0 = diag([0.3; 0.3; 0.3;   0.2; 0.2; 0.2;    1]);
    
    elseif (SizeCOM == 4 && Ts == 0.020)
    	mw0 = [-1.0; -0.5; -1.0;     -2; -1; -2;       8];
    	Sw0 = diag([0.3; 0.3; 0.3;   0.2; 0.2; 0.2;    1]);
        
    elseif (SizeCOM == 4 && Ts == 0.025)
    	mw0 = [-0.5; -0.15; -0.5;     -1.7; -0.6; -1.7;    12];
    	Sw0 = diag([0.04; 0.04; 0.04;  0.08; 0.08; 0.08;   0.5]);
        
    elseif (SizeCOM == 7 && Ts == 0.010)
        % Used for S5 10_1
        %mw0 = [-1.2 -0.98 -1.8    -3.1 -2.9 -3.2  16]';
        %Sw0 = diag([0.1; 0.1; 0.1;   0.05; 0.05; 0.05;     2]);
        % Used for the rest
        mw0 = [-1.5; -1.0; -1.5;    -3.0; -2.5; -3.0;  20];
        Sw0 = diag([0.5; 0.5; 0.5;   0.5; 0.5; 0.5;     2]);
        
    elseif (SizeCOM == 7 && Ts == 0.015)
        mw0 = [-1.5; -1.0; -1.5;    -2.5; -2.0; -2.5;  15];
        Sw0 = diag([0.2; 0.2; 0.2;   0.1; 0.1; 0.1;     1]);
    
    elseif (SizeCOM == 7 && Ts == 0.020)
        mw0 = [-1.1; -0.7; -1.0;       -2.4; -1.8; -2.4;     25];
        Sw0 = diag([0.08; 0.08; 0.08;   0.07; 0.07; 0.07;     1]);
    
    elseif (SizeCOM == 7 && Ts == 0.025)
        mw0 = [-0.60; -0.50; -0.70;     -1.7; -1.35; -1.9;   40];
        Sw0 = diag([0.02; 0.02; 0.02;   0.01; 0.01; 0.01;     1]);
        UseLambda = 0;
        
    end
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
        while any(thetai(1:6)>0) || thetai(7) < 0
            thetai = mvnrnd(mw, SwL)';
        end

        theta = struct();
        theta.stiffness = thetai(1:3)'*100;
        theta.damping = thetai(4:6)';
        theta.z_sum = thetai(7)/100;

        fprintf([' Theta: [',num2str(thetai'.*ThMask',5),']\n']);
        
        % To compare against real cloth data
        [Rwd, AllSt] = sim_ol_theta_realcloth(theta, opts);
        
        % To compare against a nonlinear model (old/new)
        %[Rwd, AllSt] = sim_ol_theta_nlmdl_new(theta, opts);
        %[Rwd, AllSt] = sim_ol_theta_nlmdl_old(theta, opts);
        
        Rwd = max(Rwd, minRwd);

        wghts_ep(:,i) = thetai;
        rwrds_ep(i) = Rwd;
        XR{1,i,epoch} = AllSt;
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
        mw = sum(weights2'.*dw)'/sum(dw);
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

%{
- Matrix of SOM states: All_StV
- Matrix of reduced SOM states: All_StVrd
- Matrix of SOM control actions: All_uSOM
- Matrix of linear control actions: All_ulin
- Matrix of simulated COM states: All_StC
- Reward: Rwd
%}


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
    
    MW2D_Thi = (MW2D(:,epoch).*ThMask)';
    
    theta = struct();
    theta.stiffness = MW2D_Thi(1:3);
    theta.damping = MW2D_Thi(4:6);
    theta.z_sum = MW2D_Thi(7);
    
    [Rwd, AllSt] = simulation_ol_theta_v3(theta, opts);
    %[Rwd, AllSt] = simulation_ol_theta_v2(theta, opts);
    %[Rwd, AllSt] = simulation_ol_theta_v1(theta, opts);
    %[Rwd, AllSt] = simulation_cl_theta(theta, opts);
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
hf1.Position = hf1.Position - [250 0 0 0];

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
hf2.Position = hf1.Position + [500 0 0 0];


full_exp_plot(dirname, 2);
if (min(min(hf2.Children.YLim)) < minRwd)
    ylim([minRwd, min(max(RWMW)+10,0)]);
end


%% 3D & 2D Plots
tSim = 0:Ts:size(AllSt.SOM,2)*Ts-Ts;
SOM_node_ctrl = [opts.nSOM*(opts.nSOM-1)+1, opts.nSOM*opts.nSOM];
SOM_coord_ctrl = [SOM_node_ctrl SOM_node_ctrl+opts.nSOM^2 SOM_node_ctrl+2*opts.nSOM^2];
COM_node_ctrl = [opts.nCOM*(opts.nCOM-1)+1, opts.nCOM*opts.nCOM];
COM_coord_ctrl = [COM_node_ctrl, COM_node_ctrl+opts.nCOM^2, COM_node_ctrl+2*opts.nCOM^2];
coord_l  = [1 opts.nCOM 1+opts.nCOM^2 opts.nCOM^2+opts.nCOM 2*opts.nCOM^2+1 2*opts.nCOM^2+opts.nCOM]; 
coord_nl = [1 opts.nSOM 1+opts.nSOM^2 opts.nSOM^2+opts.nSOM 2*opts.nSOM^2+1 2*opts.nSOM^2+opts.nSOM];


% 3D trajectories (uc & lc)
if(Plot3DTraj==1)
    hf3 = figure(3);
    hf3.Color = [1,1,1];

    plot3(AllSt.SOM(SOM_coord_ctrl(1),:), AllSt.SOM(SOM_coord_ctrl(3),:), AllSt.SOM(SOM_coord_ctrl(5),:))
    hold on
    box on; grid on;
    plot3(AllSt.SOM(SOM_coord_ctrl(2),:), AllSt.SOM(SOM_coord_ctrl(4),:), AllSt.SOM(SOM_coord_ctrl(6),:))
    plot3(AllSt.COM(COM_coord_ctrl(1),:), AllSt.COM(COM_coord_ctrl(3),:), AllSt.COM(COM_coord_ctrl(5),:), '--')
    plot3(AllSt.COM(COM_coord_ctrl(2),:), AllSt.COM(COM_coord_ctrl(4),:), AllSt.COM(COM_coord_ctrl(6),:), '--')
    plot3(AllSt.SOM(coord_nl(1),:), AllSt.SOM(coord_nl(3),:), AllSt.SOM(coord_nl(5),:))
    plot3(AllSt.SOM(coord_nl(2),:), AllSt.SOM(coord_nl(4),:), AllSt.SOM(coord_nl(6),:))
    plot3(AllSt.COM(coord_l(1),:), AllSt.COM(coord_l(3),:), AllSt.COM(coord_l(5),:), '--')
    plot3(AllSt.COM(coord_l(2),:), AllSt.COM(coord_l(4),:), AllSt.COM(coord_l(6),:), '--')
    hold off
    axis equal
end

% 2D evolutions
if(Plot2DTraj==1)
    hf4 = figure(4);
    hf4.Color = [1,1,1];

    subplot(10,2,1:2:18);
    plot(tSim, AllSt.SOM(coord_nl([1,3,5]),:));
    hold on
    plot(tSim, AllSt.COM(coord_l([1,3,5]),:), '--')
    hold off
    grid on
    set(gca,'TickLabelInterpreter','latex');
    ylim41 = ylim;

    subplot(10,2,2:2:18);
    pa4som = plot(tSim, AllSt.SOM(coord_nl([2,4,6]),:));
    hold on
    pa4com = plot(tSim, AllSt.COM(coord_l([2,4,6]),:), '--');
    hold off
    grid on
    set(gca,'TickLabelInterpreter','latex');
    ylim42 = ylim;

    ylim4 = [min(ylim41(1), ylim42(1)), max(ylim41(2), ylim42(2))];
    subplot(10,2,1:2:18);
    ylim(ylim4);
    subplot(10,2,2:2:18);
    ylim(ylim4);

    Lgnd4 = legend([pa4som(1), pa4com(1), pa4som(2), pa4com(2), pa4som(3), pa4com(3)], ...
                   '$x_{SOM}$','$x_{COM}$','$y_{SOM}$', ...
                   '$y_{COM}$','$z_{SOM}$','$z_{COM}$', ...
                   'NumColumns',3, 'Interpreter', 'latex');
    Lgnd4.Position(1) = 0.5-Lgnd4.Position(3)/2;
    Lgnd4.Position(2) = 0.04;

end


%% Update LUT for Parameters
DataTable = readtable('LearntMdl_Data.csv');
Tbl_Exp_id = (categorical(DataTable.ExpSetName) == ExpSet) & ...
             (DataTable.NExp == NExp) & (DataTable.NTrial == NTrial) & ...
             (DataTable.Ts == Ts) & (DataTable.NEpochs == NEpochs) & ...
             (DataTable.NSamples == NSamples) & (DataTable.MinRwd == minRwd) & ...
             (categorical(DataTable.Simulation) == SimType);
Tbl_Exp = DataTable(Tbl_Exp_id, :);

if (size(Tbl_Exp,1) > 1)
    error("There are multiple rows with same experiment parameters.");
elseif (size(Tbl_Exp,1) == 1)
    % Update experiment row
    Tbl_row = find(Tbl_Exp_id);
    DataTable(Tbl_row,'LastEpoch') = {e0+NEpochs};
    DataTable(Tbl_row,'Th_stiffness_x') = {ThLearnt(1)};
    DataTable(Tbl_row,'Th_stiffness_y') = {ThLearnt(2)};
    DataTable(Tbl_row,'Th_stiffness_z') = {ThLearnt(3)};
    DataTable(Tbl_row,'Th_damping_x') = {ThLearnt(4)};
    DataTable(Tbl_row,'Th_damping_y') = {ThLearnt(5)};
    DataTable(Tbl_row,'Th_damping_z') = {ThLearnt(6)};
    DataTable(Tbl_row,'Th_z_sum') = {ThLearnt(7)};
    DataTable(Tbl_row,'Rwd') = {RWMW(end)};
else
    % Add experiment row
    Tbl_row = size(DataTable,1)+1;
    DataTable(Tbl_row,'ExpSetN') = {ExpSetN};
    DataTable(Tbl_row,'ExpSetName') = {ExpSet};
    DataTable(Tbl_row,'ExpDate') = {ExpDate};
    DataTable(Tbl_row,'NExp') = {NExp};
    DataTable(Tbl_row,'NTrial') = {NTrial};
    DataTable(Tbl_row,'Ts') = {Ts};
    DataTable(Tbl_row,'nSOM') = {SizeSOM};
    DataTable(Tbl_row,'nCOM') = {SizeCOM};
    DataTable(Tbl_row,'LastEpoch') = {e0+NEpochs};
    DataTable(Tbl_row,'NEpochs') = {NEpochs};
    DataTable(Tbl_row,'NSamples') = {NSamples};
    DataTable(Tbl_row,'Simulation') = {SimType};
    DataTable(Tbl_row,'MinRwd') = {minRwd};
    
    DataTable(Tbl_row,'Th_stiffness_x') = {ThLearnt(1)};
    DataTable(Tbl_row,'Th_stiffness_y') = {ThLearnt(2)};
    DataTable(Tbl_row,'Th_stiffness_z') = {ThLearnt(3)};
    DataTable(Tbl_row,'Th_damping_x') = {ThLearnt(4)};
    DataTable(Tbl_row,'Th_damping_y') = {ThLearnt(5)};
    DataTable(Tbl_row,'Th_damping_z') = {ThLearnt(6)};
    DataTable(Tbl_row,'Th_z_sum') = {ThLearnt(7)};
    DataTable(Tbl_row,'Rwd') = {RWMW(end)};
end
DataTableCpp = DataTable(:,{'ExpSetN','NExp','NTrial','Ts','nSOM','nCOM', ...
                          'Th_stiffness_x','Th_stiffness_y','Th_stiffness_z', ...
                          'Th_damping_x','Th_damping_y','Th_damping_z', ...
                          'Th_z_sum'});
writetable(DataTable,'LearntMdl_Data.csv');
writematrix(table2array(DataTableCpp),'LearntMdl_Dcpp.csv');
fprintf('Updated LearntMdl_Data.csv and LearntMdl_Dcpp.csv\n');

