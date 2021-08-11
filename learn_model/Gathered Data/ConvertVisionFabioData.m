clear; clc;

% Trajectory
NRange = 60;
NRep = 1;


% Processing parameters
ExpSetN = 11;             % Set of experiments this data will be the base of
ExpDate = '2021_07_05';   % Date of the experiments (data folder)
SizeMesh = 8;             % Size of the mesh from the Vision node
CamTh = 0.04;             % Omit large deviations
nwght = [0.7, 0.4, 0.1];  % Weight for prev/next points (decreasing)
RegTs = 0.010;            % Regular sample time of resulting trajectory
TrimTh = 0.001;           % Trim low movements at ini/end
TrimPad = 50;             % Pad with static points
SizeSOM = 10;             % New SOM mesh side size


% Plot and save variables
animOGData = 0;
animRszData = 0;
plot3DType = 1;
plot2DType = 1;
saveMat = 0;


% Add directories
addpath(['./Session_',ExpDate,'/Converted CSV/']);
addpath(['./Session_',ExpDate,'/Input Trajectories/']);


% Calibration
calib_campose = readtable(['n',num2str(SizeMesh),'_calibration_campose.csv']);
calib_PoseRC = table2array(calib_campose(end,2:end));
calib_pRC = calib_PoseRC(1:3);
calib_qRC = [calib_PoseRC(7) calib_PoseRC(4:6)]; %wxyz
calib_RRC = quat2rotm(calib_qRC);
TRC = [calib_RRC calib_pRC'; 0 0 0 1];


% Load real trajectory data
traj_table = readtable(['range',num2str(NRange), ...
                        '_rep', num2str(NRep),'_mesh.csv']);

nSOM = traj_table.field_size(1);
if (nSOM ~= SizeMesh) 
    error("Size of the Mesh inside the data file does not match with its name");
end
TableOffset=5;
TCPOffset = -0.03; % IK on WAM computed with d7=6cm, it's actually 9cm
Npt = size(traj_table,1);
t0 = table2array(traj_table(1,3))/1e9;


% Load real TCP trajectory data if available
try
    tcp_table = readtable(['range',num2str(NRange), ...
                           '_rep', num2str(NRep),'_tcp.csv']);
    tcp_t0 = table2array(tcp_table(1,3))/1e9;
    tcp_data = table2array(tcp_table(:, [3,5:end]));
    tcp_data(:,1) = tcp_data(:,1)/1e9 - t0;
    
    use_tcp = 1;
catch MErr
    use_tcp = 0;
end

% Load converted SOMstate data if possible
try
    somstate_table = readtable(['range',num2str(NRange), ...
                                '_rep', num2str(NRep),'_somstate.csv']);
    somstate_t0 = table2array(somstate_table(1,3))/1e9;
    somstate_data = table2array(somstate_table(:, [3,6:end]));
    somstate_data(:,1) = somstate_data(:,1)/1e9 - t0;
    
    use_somstate = 1;
catch MErr
    use_somstate = 0;
end

% Load input trajectory
InputTraj_table = readtable(['trajectory_range_',num2str(NRange), ...
                             '_rep_', num2str(NRep),'.csv']);
InputTraj = table2array(InputTraj_table(:,9:15)); %xyzq, q=wxyz
InputTrajT = table2array(InputTraj_table(:,1));
InputTrajTs = InputTrajT(2)-InputTrajT(1);


% Corner coordinates
coord_nl = [1 nSOM 1+nSOM^2 nSOM+nSOM^2 1+2*nSOM^2 nSOM+2*nSOM^2];
coord_ctrl = [nSOM^2-nSOM+1 nSOM^2 2*nSOM^2-nSOM+1 2*nSOM^2 3*nSOM^2-nSOM+1 3*nSOM^2];


% Plot mesh
animAny = animOGData || animRszData;
if (animAny)
    fig1 = figure(1);
    fig1.Color = [1,1,1];
end

lt = [-0.5 0.6 -1 0.1 -0.3 0.6];
pov = [-5 10]; %[-25 20]; %[-70 20];

All_SOMpos = zeros(3*nSOM^2, Npt);
All_SOMt = zeros(1, Npt);
for ptk=1:Npt
    tCAM = table2array(traj_table(ptk, 3))/1e9 - t0;
    xCAM = table2array(traj_table(ptk, TableOffset+1:TableOffset+nSOM^2))';
    yCAM = table2array(traj_table(ptk, TableOffset+nSOM^2+1:TableOffset+2*nSOM^2))';
    zCAM = table2array(traj_table(ptk, TableOffset+2*nSOM^2+1:TableOffset+3*nSOM^2))';

    dxC = max(xCAM)-min(xCAM);
    dyC = max(yCAM)-min(yCAM);

    xCsq = reshape(xCAM, nSOM,nSOM);
    yCsq = reshape(yCAM, nSOM,nSOM);
    zCsq = reshape(zCAM, nSOM,nSOM);

    % Desired: x increasing, then y decreasing
    if range(xCsq(:,1)) < 0.4*dxC
        % First 10 on same x --> Transpose!
        xCsq = xCsq';
        yCsq = yCsq';
        zCsq = zCsq';
    end
    if (xCsq(end,1) - xCsq(1,1)) < 0
        % Decreasing x --> Invert!
        xCsq = xCsq(end:-1:1,:);
        yCsq = yCsq(end:-1:1,:);
        zCsq = zCsq(end:-1:1,:);
    end
    if (yCsq(1,end) - yCsq(1,1)) > 0
        % Increasing y --> Invert!
        xCsq = xCsq(:,end:-1:1);
        yCsq = yCsq(:,end:-1:1);
        zCsq = zCsq(:,end:-1:1);
    end

    xCAMord = reshape(xCsq, nSOM^2,1);
    yCAMord = reshape(yCsq, nSOM^2,1);
    zCAMord = reshape(zCsq, nSOM^2,1);

    % Apply Transform to all points
    xSOM = zeros(size(xCAMord));
    ySOM = zeros(size(yCAMord));
    zSOM = zeros(size(zCAMord));
    for pj=1:length(xCAMord)

        % Filter NaNs in depth (z)
        if isnan(zCAMord(pj))
            zSum = 0; zN = 0; pj1 = pj-1; 
            if (floor(pj1/nSOM)>0) && ~isnan(zCAMord(pj-nSOM))
            	zSum=zSum+zCAMord(pj-nSOM); zN=zN+1;
            end
            if(floor(pj1/nSOM)<nSOM-1) && ~isnan(zCAMord(pj+nSOM))
                zSum=zSum+zCAMord(pj+nSOM); zN=zN+1;
            end
            if(mod(pj1,nSOM)>0) && ~isnan(zCAMord(pj-1))
                zSum=zSum+zCAMord(pj-1); zN=zN+1;
            end
            if(mod(pj1,nSOM)<nSOM-1) && ~isnan(zCAMord(pj+1))
                zSum=zSum+zCAMord(pj+1); zN=zN+1;
            end

            zCAMord(pj) = zSum/zN;
        end

        Pt = [xCAMord(pj) yCAMord(pj) zCAMord(pj) 1]';

        RPt = TRC*Pt;

        xSOM(pj) = RPt(1);
        ySOM(pj) = RPt(2);
        zSOM(pj) = RPt(3);
    end
    SOMpos = [xSOM;ySOM;zSOM];
    
    % Store data
    All_SOMpos(:,ptk) = SOMpos;
    All_SOMt(:,ptk) = tCAM;

    % Plot Data in SOM coordinates
    if (animOGData > 0)
        scatter3(xSOM(1:nSOM), ySOM(1:nSOM), zSOM(1:nSOM), '.r');
        hold on
        scatter3(xSOM(1), ySOM(1), zSOM(1), '*r');
        scatter3(xSOM(nSOM+1:end), ySOM(nSOM+1:end), zSOM(nSOM+1:end), '.b');
        hold off
        legend off
        axis equal;
        axis(lt);
        fig1.Children.View = pov;

        box on;
        xlabel('X', 'Interpreter','latex');
        ylabel('Y', 'Interpreter','latex');
        zlabel('Z', 'Interpreter','latex');

        set(gca, 'TickLabelInterpreter','latex');
        pause(1e-6);
    end

end


%% Omit large deviations
dSOMpos = diff(All_SOMpos,1,2);
dSOMxyz = reshape(permute(dSOMpos, [1,3,2]), [nSOM^2,3,Npt-1]);
dSOMnorm = permute(vecnorm(dSOMxyz,2,2), [1,3,2]);

dSOMth = dSOMnorm >= CamTh;
dSOMth = [zeros(nSOM^2,1), dSOMth];
dSOMth3 = [dSOMth;dSOMth;dSOMth];

SOMpos_aux = ([zeros(3*nSOM^2,1) All_SOMpos(:,1:end-1).*dSOMth3(:,2:end)] + ...
              [All_SOMpos(:,2:end).*dSOMth3(:,1:end-1) zeros(3*nSOM^2,1)])/2 + ...
              All_SOMpos.*~dSOMth3;
          
dSOMpos_aux = diff(SOMpos_aux,1,2);
dSOMxyz_aux = reshape(permute(dSOMpos_aux,[1,3,2]), [nSOM^2,3,Npt-1]);
dSOMnorm_aux = permute(vecnorm(dSOMxyz_aux,2,2),[1,3,2]);
dSOMth_aux = dSOMnorm_aux>=CamTh;


%% Filter with neighbor means
fwidth = size(nwght,2);
ffull = [nwght(end:-1:1), 1, nwght];
% Choose one
%{
SOMpos_ext = [All_SOMpos(:,1)*ones(1,fwidth), All_SOMpos, ...
              All_SOMpos(:,end)*ones(1,fwidth)];
%}
%{x
SOMpos_ext = [SOMpos_aux(:,1)*ones(1,fwidth), SOMpos_aux, ...
              SOMpos_aux(:,end)*ones(1,fwidth)];
%}
SOMpos_f = zeros(size(All_SOMpos));
for fwi=1:length(ffull)
    SOMpos_fi = ffull(fwi)*SOMpos_ext(:,fwi:end-length(ffull)+fwi);
    SOMpos_f = SOMpos_f + SOMpos_fi;
end
SOMpos_f = SOMpos_f/sum(ffull);

%SOMpos_f = (nwght*SOMpos_ext(:,1:end-2) + nwght*SOMpos_ext(:,3:end) + ...
%             SOMpos_ext(:,2:end-1))/(1+2*nwght);


%% Interpolate data for a fixed Ts
SOMt_end = All_SOMt(end);
Reg_SOMt = 0:RegTs:SOMt_end;
Reg_Npt = length(Reg_SOMt);
Reg_SOMpos = zeros(3*nSOM^2, Reg_Npt);

% Select one data series to regularize
%SOMpos = All_SOMpos;
%SOMpos = SOMpos_aux;
SOMpos = SOMpos_f;

for Rpt = 1:Reg_Npt
    Reg_ti = Reg_SOMt(Rpt);
    SOM_ti = sum(All_SOMt <= Reg_ti);
    SOMt_prev = All_SOMt(SOM_ti);
    SOMt_next = All_SOMt(SOM_ti+1);
    if SOMt_prev < SOMt_end
        t_lerp = (Reg_ti - SOMt_prev)/(SOMt_next - SOMt_prev);
        SOMposRi = t_lerp*(SOMpos(:,SOM_ti+1) - SOMpos(:,SOM_ti)) ...
                  + SOMpos(:,SOM_ti);
        Reg_SOMpos(:,Rpt) = SOMposRi;
    else
        % This shouldn't happen by construction, but kept as a failsafe
        t_lerp = 0;
        SOMposRi = SOMpos(:,SOM_ti);
        Reg_SOMpos(:,Rpt) = SOMposRi;
    end

end


%% Interpolate to new SOM size
if SizeSOM == nSOM
    Rsz_SOMpos = Reg_SOMpos;
    Rsz_SOMt = Reg_SOMt;
    Rsz_coord_nl = coord_nl;
    Rsz_coord_ctrl = coord_ctrl;
else
    Reg_SOMx = Reg_SOMpos(1:nSOM^2,:);
    Reg_SOMy = Reg_SOMpos(nSOM^2+1:2*nSOM^2,:);
    Reg_SOMz = Reg_SOMpos(2*nSOM^2+1:3*nSOM^2,:);
    Rsz_SOMx = zeros(SizeSOM^2, Reg_Npt);
    Rsz_SOMy = zeros(SizeSOM^2, Reg_Npt);
    Rsz_SOMz = zeros(SizeSOM^2, Reg_Npt);

    Reg_SOMlsp = 0:1/(nSOM-1):1;
    Rsz_SOMlsp = 0:1/(SizeSOM-1):1;

    for pti = 1:SizeSOM^2
        ptiRow = floor((pti-1)/SizeSOM) + 1;
        ptiCol = mod(pti-1, SizeSOM) + 1;

        Reg_prevRow = sum(Reg_SOMlsp <= Rsz_SOMlsp(ptiRow));
        Reg_prevCol = sum(Reg_SOMlsp <= Rsz_SOMlsp(ptiCol));

        % (pt-prev)/(next-prev) = (pt-prev)/step = (pt-prev)*spaces
        lerpRow = (Rsz_SOMlsp(ptiRow) - Reg_SOMlsp(Reg_prevRow))*(nSOM-1);
        lerpCol = (Rsz_SOMlsp(ptiCol) - Reg_SOMlsp(Reg_prevCol))*(nSOM-1);

        % 2D interpolation: i, i+1, i+n, i+n+1
        idx_min = (Reg_prevRow-1)*nSOM + Reg_prevCol;

        Rsz_SOMx(pti,:) = Reg_SOMx(idx_min, :)*(1-lerpRow)*(1-lerpCol);
        Rsz_SOMy(pti,:) = Reg_SOMy(idx_min, :)*(1-lerpRow)*(1-lerpCol);
        Rsz_SOMz(pti,:) = Reg_SOMz(idx_min, :)*(1-lerpRow)*(1-lerpCol);
        if (lerpCol > 0)
            Rsz_SOMx(pti,:) = Rsz_SOMx(pti,:) + Reg_SOMx(idx_min+1, :)*(1-lerpRow)*(lerpCol);
            Rsz_SOMy(pti,:) = Rsz_SOMy(pti,:) + Reg_SOMy(idx_min+1, :)*(1-lerpRow)*(lerpCol);
            Rsz_SOMz(pti,:) = Rsz_SOMz(pti,:) + Reg_SOMz(idx_min+1, :)*(1-lerpRow)*(lerpCol);
        end
        if (lerpRow > 0)
            Rsz_SOMx(pti,:) = Rsz_SOMx(pti,:) + Reg_SOMx(idx_min+nSOM, :)*(lerpRow)*(1-lerpCol);
            Rsz_SOMy(pti,:) = Rsz_SOMy(pti,:) + Reg_SOMy(idx_min+nSOM, :)*(lerpRow)*(1-lerpCol);
            Rsz_SOMz(pti,:) = Rsz_SOMz(pti,:) + Reg_SOMz(idx_min+nSOM, :)*(lerpRow)*(1-lerpCol);
        end
        if (lerpRow > 0 && lerpCol > 0)
            Rsz_SOMx(pti,:) = Rsz_SOMx(pti,:) + Reg_SOMx(idx_min+nSOM+1, :)*(lerpRow)*(lerpCol);
            Rsz_SOMy(pti,:) = Rsz_SOMy(pti,:) + Reg_SOMy(idx_min+nSOM+1, :)*(lerpRow)*(lerpCol);
            Rsz_SOMz(pti,:) = Rsz_SOMz(pti,:) + Reg_SOMz(idx_min+nSOM+1, :)*(lerpRow)*(lerpCol);
        end

    end

    Rsz_SOMpos = [Rsz_SOMx; Rsz_SOMy; Rsz_SOMz];
    Rsz_SOMt = Reg_SOMt;
    Rsz_coord_nl = [1             SizeSOM, ...
                    1+1*SizeSOM^2 SizeSOM+1*SizeSOM^2, ...
                    1+2*SizeSOM^2 SizeSOM+2*SizeSOM^2];
    Rsz_coord_ctrl = [1*SizeSOM^2-SizeSOM+1 1*SizeSOM^2, ...
                      2*SizeSOM^2-SizeSOM+1 2*SizeSOM^2, ...
                      3*SizeSOM^2-SizeSOM+1 3*SizeSOM^2];

    if (animRszData > 0)
        for Rpt = 1:Reg_Npt
            SOMposRi = Rsz_SOMpos(:,Rpt);

            Rsz_xSOMi = SOMposRi(1:SizeSOM^2);
            Rsz_ySOMi = SOMposRi(SizeSOM^2+1:2*SizeSOM^2);
            Rsz_zSOMi = SOMposRi(2*SizeSOM^2+1:3*SizeSOM^2);

            scatter3(Rsz_xSOMi(1:SizeSOM), Rsz_ySOMi(1:SizeSOM), Rsz_zSOMi(1:SizeSOM), '.r');
            hold on
            scatter3(Rsz_xSOMi(1), Rsz_ySOMi(1), Rsz_zSOMi(1), '*r');
            scatter3(Rsz_xSOMi(SizeSOM+1:end), Rsz_ySOMi(SizeSOM+1:end), Rsz_zSOMi(SizeSOM+1:end), '.b');
            hold off
            legend off
            axis equal;
            axis(lt);
            fig1.Children.View = pov;

            box on;
            xlabel('X', 'Interpreter','latex');
            ylabel('Y', 'Interpreter','latex');
            zlabel('Z', 'Interpreter','latex');

            set(gca, 'TickLabelInterpreter','latex');
            pause(1e-6);

        end
    end
end


%% Trim start and end
Rsz_dSOMpos = diff(Rsz_SOMpos,1,2);
Rsz_dSOMxyz = reshape(permute(Rsz_dSOMpos, [1,3,2]), [SizeSOM^2,3,Reg_Npt-1]);
Rsz_dSOMnorm = permute(vecnorm(Rsz_dSOMxyz,2,2), [1,3,2]);
Rsz_dSOMth = mean(Rsz_dSOMnorm)' > TrimTh;
Rsz_SOMlims = [find(Rsz_dSOMth,1,'first'), find(Rsz_dSOMth,1,'last')];

Proc_SOMpos = [Rsz_SOMpos(:,Rsz_SOMlims(1)).*ones(1,TrimPad), ...
               Rsz_SOMpos(:,Rsz_SOMlims(1):Rsz_SOMlims(2)), ...
               Rsz_SOMpos(:,Rsz_SOMlims(2)).*ones(1,TrimPad)];

Proc_SOMtx = [Rsz_SOMt(:,Rsz_SOMlims(1))-TrimPad*RegTs:RegTs:Rsz_SOMt(:,Rsz_SOMlims(1))-RegTs, ...
              Rsz_SOMt(:,Rsz_SOMlims(1):Rsz_SOMlims(2)), ...
              Rsz_SOMt(:,Rsz_SOMlims(2))+RegTs:RegTs:Rsz_SOMt(:,Rsz_SOMlims(2))+TrimPad*RegTs];
Proc_SOMt = Proc_SOMtx - min(Proc_SOMtx);


%% Add velocity
Proc_SOMvel = (Proc_SOMpos(:,2:end) - Proc_SOMpos(:,1:end-1))/RegTs;
Proc_SOMvel = [zeros(3*SizeSOM^2,1) Proc_SOMvel];
Proc_SOMst  = [Proc_SOMpos; Proc_SOMvel];


%% Plot trajectories

% In 3D space
if (plot3DType > 0)
    fig2 = figure(2);
    fig2.Color = [1,1,1];
    plot3(InputTraj(:,1), InputTraj(:,2), InputTraj(:,3), '-k')
    hold on
    
    plot3(All_SOMpos(coord_nl(1:2),:)',All_SOMpos(coord_nl(3:4),:)', ...
          All_SOMpos(coord_nl(5:6),:)',':b')
    plot3(All_SOMpos(coord_ctrl(1:2),:)',All_SOMpos(coord_ctrl(3:4),:)', ...
          All_SOMpos(coord_ctrl(5:6),:)',':b')
          
    plot3(Rsz_SOMpos(Rsz_coord_nl(1:2),:)',Rsz_SOMpos(Rsz_coord_nl(3:4),:)', ...
          Rsz_SOMpos(Rsz_coord_nl(5:6),:)','-','Color',[1 0.5 0])
    plot3(Rsz_SOMpos(Rsz_coord_ctrl(1:2),:)',Rsz_SOMpos(Rsz_coord_ctrl(3:4),:)', ...
          Rsz_SOMpos(Rsz_coord_ctrl(5:6),:)','-','Color',[1 0.5 0])
      
    plot3(Proc_SOMpos(Rsz_coord_nl(1:2),:)',Proc_SOMpos(Rsz_coord_nl(3:4),:)', ...
          Proc_SOMpos(Rsz_coord_nl(5:6),:)','-','Color',[0 0.6 0])
    plot3(Proc_SOMpos(Rsz_coord_ctrl(1:2),:)',Proc_SOMpos(Rsz_coord_ctrl(3:4),:)', ...
          Proc_SOMpos(Rsz_coord_ctrl(5:6),:)','-','Color',[0 0.6 0])
    
    if (use_tcp == 1)
        plot3(tcp_data(:,2), tcp_data(:,3), tcp_data(:,4), '--m')
    end
      
    hold off
    axis equal; grid on; box on
    set(gca,'TickLabelInterpreter','latex');
    
end



% Evolution of each coordinate over time
if (plot2DType > 0)
    fig3 = figure(3);
    fig3.Color = [1,1,1];

    subplot(2,2,1)
    plot(All_SOMt, All_SOMpos(coord_nl([1,3,5]),:)',':b','LineWidth',1)
    hold on
    plot(Rsz_SOMt, Rsz_SOMpos(Rsz_coord_nl([1,3,5]),:)','-','Color',[1 0.5 0],'LineWidth',1)
    plot(Proc_SOMtx, Proc_SOMpos(Rsz_coord_nl([1,3,5]),:)','-','Color',[0 0.6 0],'LineWidth',1)
    hold off
    grid on
    set(gca,'TickLabelInterpreter','latex');
    
    subplot(2,2,2)
    plot(All_SOMt, All_SOMpos(coord_nl([2,4,6]),:)',':b','LineWidth',1)
    hold on
    plot(Rsz_SOMt, Rsz_SOMpos(Rsz_coord_nl([2,4,6]),:)','-','Color',[1 0.5 0],'LineWidth',1)
    plot(Proc_SOMtx, Proc_SOMpos(Rsz_coord_nl([2,4,6]),:)','-','Color',[0 0.6 0],'LineWidth',1)
    hold off
    grid on
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,3)
    plot(All_SOMt, All_SOMpos(coord_ctrl([1,3,5]),:)',':b','LineWidth',1)
    hold on
    plot(Rsz_SOMt, Rsz_SOMpos(Rsz_coord_ctrl([1,3,5]),:)','-','Color',[1 0.5 0],'LineWidth',1)
    plot(Proc_SOMtx, Proc_SOMpos(Rsz_coord_ctrl([1,3,5]),:)','-','Color',[0 0.6 0],'LineWidth',1)
    hold off
    grid on
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,4)
    plot(All_SOMt, All_SOMpos(coord_ctrl([2,4,6]),:)',':b','LineWidth',1)
    hold on
    plot(Rsz_SOMt, Rsz_SOMpos(Rsz_coord_ctrl([2,4,6]),:)','-','Color',[1 0.5 0],'LineWidth',1)
    plot(Proc_SOMtx, Proc_SOMpos(Rsz_coord_ctrl([2,4,6]),:)','-','Color',[0 0.6 0],'LineWidth',1)
    hold off
    grid on
    set(gca,'TickLabelInterpreter','latex');

end


%% Somstate plots
if (use_somstate == 1)
    fig4 = figure(4);
    fig4.Color = [1,1,1];
    
    somstate_ctrl_data = somstate_data(:,[1,coord_ctrl+1]);
    somstate_nl_data = somstate_data(:,[1,coord_nl+1]);

    plot3(somstate_ctrl_data(:,2:3),somstate_ctrl_data(:,4:5),somstate_ctrl_data(:,6:7))
    grid on; box on; axis equal
    hold on
    plot3(somstate_nl_data(:,2:3),somstate_nl_data(:,4:5),somstate_nl_data(:,6:7))
    hold off
    set(gca, 'TickLabelInterpreter','latex');
    xlabel('X', 'Interpreter','latex');
    ylabel('Y', 'Interpreter','latex');
    zlabel('Z', 'Interpreter','latex');

    hold on
    plot3(tcp_data(:,2), tcp_data(:,3), tcp_data(:,4), '--m')
    hold off

    
    fig5 = figure(5);
    fig5.Color = [1,1,1];

    subplot(2,2,1)
    plot(somstate_ctrl_data(:,1), somstate_ctrl_data(:,[2,4,6]),'LineWidth',1)
    grid on
    set(gca,'TickLabelInterpreter','latex');
    
    subplot(2,2,2)
    plot(somstate_ctrl_data(:,1), somstate_ctrl_data(:,[3,5,7]),'LineWidth',1)
    grid on
    set(gca,'TickLabelInterpreter','latex');
    
    subplot(2,2,3)
    plot(somstate_nl_data(:,1), somstate_nl_data(:,[2,4,6]),'LineWidth',1)
    grid on
    set(gca,'TickLabelInterpreter','latex');
    
    subplot(2,2,4)
    plot(somstate_nl_data(:,1), somstate_nl_data(:,[3,5,7]),'LineWidth',1)
    grid on
    set(gca,'TickLabelInterpreter','latex');
end



%% Save new data
SOMrec = struct();
SOMrec.states = Proc_SOMst;
SOMrec.time = Proc_SOMt;

dirname = ['Session_',ExpDate,'/Processed MAT S', num2str(ExpSetN)];
newfilename = ['SOMtraj',num2str(NRange),'_',num2str(NRep),'_',num2str(RegTs*1000),'ms.mat'];

if (~isfolder(dirname) && saveMat>0)
    mkdir(dirname);
end

if (saveMat > 0)
    save([dirname,'/',newfilename],'SOMrec');
end



