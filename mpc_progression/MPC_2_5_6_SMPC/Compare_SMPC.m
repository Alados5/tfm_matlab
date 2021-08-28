clc; clear;

%% Data

% Results for Hp=15; sigmaD=0.003; sigmaN=0.001; W_Q=0.01; W_R=1;
DMPCe = [14.0759, 14.0064, 14.8566, 14.7680, 15.0105, ...
         13.6003, 13.8963, 13.1647, 14.6818, 15.0395];
DMPCt = [09.6706, 09.9934, 09.7034, 09.9613, 09.6760, ...
         09.8435, 09.3321, 09.8466, 09.6159, 09.9736];

SMPCe = [12.8655, 13.5315, 13.7562, 13.5828, 13.6239, ...
         14.8037, 14.3465, 13.4941, 13.7849, 14.0830];
SMPCt = [09.9002, 10.1033, 10.1122, 09.9503, 09.4611, ...
         09.5086, 09.7064, 10.2092, 09.7823, 09.7820];

     
%% Compute 

ND = length(DMPCe);
NS = length(SMPCe);

muDe = mean(DMPCe);
muDt = mean(DMPCt);
sigmaDe = std(DMPCe);
sigmaDt = std(DMPCt);

muSe = mean(SMPCe);
muSt = mean(SMPCt);
sigmaSe = std(SMPCe);
sigmaSt = std(SMPCt);

Ze = abs(muSe - muDe)/sqrt((sigmaSe/sqrt(NS))^2 + (sigmaDe/sqrt(ND))^2);
Zt = abs(muSt - muDt)/sqrt((sigmaSt/sqrt(NS))^2 + (sigmaDt/sqrt(ND))^2);

fprintf(['Z-test results \n', ...
         ' - For errors: Z=',num2str(Ze), '\n', ...
         ' - For times: Z=',num2str(Zt), '\n']);

%[muDe sigmaDe muDt sigmaDt; muSe sigmaSe muSt sigmaSt]
     
xe = min(muDe,muSe)-5*max(sigmaDe,sigmaSe):.01:max(muDe,muSe)+5*max(sigmaDe,sigmaSe);
xt = min(muDt,muSt)-5*max(sigmaDt,sigmaSt):.01:max(muDt,muSt)+5*max(sigmaDt,sigmaSt);

yDe = normpdf(xe,muDe,sigmaDe);
ySe = normpdf(xe,muSe,sigmaSe);
yDt = normpdf(xt,muDt,sigmaDt);
ySt = normpdf(xt,muSt,sigmaSt);


%% Plot

fig1 = figure(1);
fig1.Color = [1,1,1];
fig1.Units = 'normalized';
fig1.Position = [0.5 0.6 0.45 0.25];

subplot(10,2,1:2:18)
plot(xe,[yDe; ySe])
grid on
set(gca, 'TickLabelInterpreter', 'latex');
xlabel('RMSE [mm]', 'Interpreter', 'latex')
%ylabel('Position [m]', 'Interpreter', 'latex')

subplot(10,2,2:2:18)
plot(xt,[yDt; ySt])
grid on
set(gca, 'TickLabelInterpreter', 'latex');
xlabel('Time [ms]', 'Interpreter', 'latex')
%ylabel('Position [m]', 'Interpreter', 'latex')

Lgnd = legend('DMPC','SMPC', 'Orientation','horizontal', ...
              'Interpreter', 'latex');
Lgnd.Position(1) = 0.5-Lgnd.Position(3)/2;
Lgnd.Position(2) = 0.01;



