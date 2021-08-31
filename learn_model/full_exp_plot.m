function full_exp_plot(dirname, hfn, plotType)

if nargin < 2
    hfn = 2;
    plotType = 1;
elseif nargin < 3
    plotType = 1;
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
    
    MWfile = load([dirname,'/MW_',epochrange,'.mat']);
    MWfile = MWfile.MW;
    MW(:,:,NEpochs*(file-1)+1:NEpochs*file+1) = MWfile;
    
    THfile = load([dirname,'/TH_',epochrange,'.mat']);
    THfile = THfile.TH;
    TH(:,:,NEpochs*(file-1)+1:NEpochs*file) = THfile;
    
    SWfile = load([dirname,'/SW_',epochrange,'.mat']);
    SWfile = SWfile.SW;
    SW(:,:,NEpochs*(file-1)+1:NEpochs*file+1) = SWfile;
end

if(plotType == 1)
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
    
elseif(plotType == 2)
    MW = permute(MW,[1,3,2])';
    
    subplot(3,1,1)
    plot(0:NFiles*NEpochs, 100*MW(:,1:3))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$k$ [N/m]', 'Interpreter','latex');
    subplot(3,1,2)
    plot(0:NFiles*NEpochs, MW(:,4:6))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$b$ [Ns/m]', 'Interpreter','latex');
    subplot(3,1,3)
    plot(0:NFiles*NEpochs, 0.01*MW(:,7))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$\Delta l_{0z}$ [m]', 'Interpreter','latex');
    
elseif(plotType == 3)
    SWD = zeros(size(SW,1),1,size(SW,3));
    for i=1:301, SWD(:,:,i) = diag(SW(:,:,i)); end
    SWD = permute(SWD,[1,3,2])';
    
    subplot(3,1,1)
    plot(0:NFiles*NEpochs, 100*SWD(:,1:3))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$k$ [N/m]', 'Interpreter','latex');
    subplot(3,1,2)
    plot(0:NFiles*NEpochs, SWD(:,4:6))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$b$ [Ns/m]', 'Interpreter','latex');
    subplot(3,1,3)
    plot(0:NFiles*NEpochs, 0.01*SWD(:,7))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$\Delta l_{0z}$ [m]', 'Interpreter','latex');
    
elseif(plotType == 4)
    MW = permute(MW,[1,3,2])';
    SWD = zeros(size(SW,1),1,size(SW,3));
    for i=1:301, SWD(:,:,i) = diag(SW(:,:,i)); end
    SWD = sqrt(permute(SWD,[1,3,2])');
    
    hf_full.Units = 'normalized';
    hf_full.Position = [0.05 0.35 0.9 0.45];
    
    subplot(2,3,1)
    plot(0:NFiles*NEpochs, 100*MW(:,1:3))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$\mu(k)$ [N/m]', 'Interpreter','latex');
    title('\textbf{Stiffness}', 'Interpreter','latex');
    subplot(2,3,2)
    plot(0:NFiles*NEpochs, MW(:,4:6))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$\mu(b)$ [Ns/m]', 'Interpreter','latex');
    title('\textbf{Damping}', 'Interpreter','latex');
    subplot(2,3,3)
    plot(0:NFiles*NEpochs, 0.01*MW(:,7))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$\mu(\Delta l_{0z})$ [m]', 'Interpreter','latex');
    title('\textbf{Length}', 'Interpreter','latex');

    subplot(2,3,4)
    plot(0:NFiles*NEpochs, 100*SWD(:,1:3))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$\sigma(k)$ [N/m]', 'Interpreter','latex');
    subplot(2,3,5)
    plot(0:NFiles*NEpochs, SWD(:,4:6))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$\sigma(b)$ [Ns/m]', 'Interpreter','latex');
    subplot(2,3,6)
    plot(0:NFiles*NEpochs, 0.01*SWD(:,7))
    grid on
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Epoch', 'Interpreter','latex');
    ylabel('$\sigma(\Delta l_{0z})$ [m]', 'Interpreter','latex');

elseif(plotType == 5)
    i = round(2/3*size(TH,3));
    rwi = RW(:,:,i);
    rwi = (rwi-min(rwi))/range(rwi);
    
    rwi(rwi==0)= min(rwi(rwi~=0));
    rwi = (rwi-min(rwi))/range(rwi);
    
    nth = size(TH,1);
    nsm = size(TH,2);

    colormap jet;
    for th=1:nth
        subplot(1,nth,th)
        scatter(i*ones(nsm,1), TH(th,:,i)', [], 1-rwi', 'filled')
        hold on
        scatter(i, MW(th,:,i)', 100, 'om', 'LineWidth',2);
        scatter(i+1, MW(th,:,i+1)', 100, 'om', 'LineWidth',2);
        hold off
        grid on
        set(gca,'TickLabelInterpreter','latex');
        xlim([i-0.5,i+1.5]);
    end


    
    
    
end

end