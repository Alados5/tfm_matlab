function full_exp_plot(dirname, hfn)

if nargin < 2
    hfn = 2;
end

hf_full = figure(hfn);
hf_full.Color = [1,1,1];

ExpFiles = dir(fullfile(['./',dirname],'RW_*.mat'));
NFiles = size(ExpFiles,1);

CellNums = split(ExpFiles(1).name(4:end-4),"-",2);
NEpochs = str2double(CellNums{2}) - str2double(CellNums{1});

for file=1:NFiles
    epochrange = [num2str(NEpochs*(file-1)), '-', num2str(NEpochs*file)];
    
    RWfile = load([dirname,'/RW_',epochrange,'.mat']);
    RWfile = RWfile.RW;
    RW(:,:,NEpochs*(file-1)+1:NEpochs*file) = RWfile;
    
    RWMWfile = load([dirname,'/RWMW_',epochrange,'.mat']);
    RWMWfile = RWMWfile.RWMW;
    RWMW(1,NEpochs*(file-1)+1:NEpochs*file+1) = RWMWfile;
end

% Evolution of sample mean rewards
RW2D = permute(RW, [2,3,1]);
RWm = mean(RW2D);

subplot(2,1,1);
plot(0:NFiles*NEpochs-1, RWm, 'LineWidth',1);
grid on
set(gca,'TickLabelInterpreter','latex');
xlim([0 NFiles*NEpochs]);
xlabel('Epoch', 'Interpreter','latex');
ylabel('Reward', 'Interpreter','latex');
title('\textbf{Mean of sample rewards}', 'Interpreter','latex');

% Evolution of resulting means

subplot(2,1,2);
plot(0:NFiles*NEpochs, RWMW, 'LineWidth',1);
hold on
plot([0 NFiles*NEpochs], [RWMW(1) RWMW(1)], '--k');
hold off
grid on
set(gca,'TickLabelInterpreter','latex');
xlabel('Epoch', 'Interpreter','latex');
ylabel('Reward', 'Interpreter','latex');
title('\textbf{Reward of resulting mean}', 'Interpreter','latex');

end