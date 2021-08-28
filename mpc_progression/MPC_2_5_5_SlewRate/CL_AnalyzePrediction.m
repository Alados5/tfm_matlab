clear; clc; close all;

%% Initialize

NTraj = 6;
Ts = 0.01;
W_Q = 0.2; %0.2 %0.05
W_R = 1;
Du = 1;
LdData = 1;

AllHps = (5:1:30)'; %(1:1:50)';
AllTTimes = zeros(length(AllHps),1);
AllMTimes = zeros(length(AllHps),1);
AllErrors = zeros(length(AllHps),1);

filename = ['hp_analysis_results\AllResults_Traj',num2str(NTraj), ...
            '_Ts',num2str(Ts*1000), '_Du',num2str(Du), '_Q',num2str(W_Q), ...
            '_R',num2str(W_R), '.mat'];


%% Main loop

if (LdData == 0)
    for i=1:length(AllHps)

        Hpi = AllHps(i);
        results = CL_SimulationFcn(NTraj, Ts, Hpi, Du, W_Q, W_R);

        AllTTimes(i) = results.TotalTime;
        AllMTimes(i) = results.MeanTime;
        AllErrors(i) = results.RMSE;

    end
    
    AllResults = struct();
    AllResults.AllTTimes = AllTTimes;
    AllResults.AllMTimes = AllMTimes;
    AllResults.AllErrors = AllErrors;
    AllResults.AllHps    = AllHps;

    save(filename,'AllResults');
end


%% Load
if (LdData==1)
    load(filename);
    AllTTimes = AllResults.AllTTimes;
    AllMTimes = AllResults.AllMTimes;
    AllErrors = AllResults.AllErrors;
    AllHps    = AllResults.AllHps;
end


%% Process

LowErr = find(AllErrors<=min(AllErrors)+0.1*range(AllErrors)) + AllHps(1)-1;
HighCT = find(AllMTimes>Ts) + AllHps(1)-1;


%% Plot

fig1 = figure(1);
fig1.Units = 'normalized';
fig1.Position = [0.475 0.6 0.5 0.3];
fig1.Color = [1,1,1];

yyaxis left
plot(AllHps, AllErrors*1000)
hold on
%plot(xlim, 1000*[AllErrors(LowErr(1)),AllErrors(LowErr(1))], '--b');
plot([LowErr(1), LowErr(1)], ylim, '--b');
hold off
ylabel('RMSE [mm]', 'Interpreter', 'latex')

yyaxis right
plot(AllHps,AllMTimes*1000)
hold on
%plot(xlim, 1000*[Ts,Ts],'--r')
plot([HighCT(1), HighCT(1)], ylim, '--r');
hold off
ylabel('Mean step time [ms]', 'Interpreter', 'latex')

grid on
set(gca, 'TickLabelInterpreter', 'latex');
xlabel('Prediction Horizon, $H_p$', 'Interpreter', 'latex')
if (Du==0)
    use_uDu = '$u$';
else
    use_uDu = '$\Delta u$';
end
title(['\textbf{Error and Computational time vs. Prediction Horizon (using ',use_uDu,')}'], 'Interpreter', 'latex')



