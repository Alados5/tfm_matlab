clear; clc;

% Parameters
Ts = 0.01;
nSOM = 10;
lCloth = 0.3;
cCloth = [0, -0.4, 0.2];
aCloth = 0;
NT = 50;
savefiles =0;

% Experiment: FROM 8 TO X
NExp = 10;

% Get initial positions (cloth assumed in XZ)
u_ini = [cCloth(1)-lCloth/2*cos(aCloth); cCloth(1)+lCloth/2*cos(aCloth);
         cCloth(2)-lCloth/2*sin(aCloth); cCloth(2)+lCloth/2*sin(aCloth);
         cCloth(3)+lCloth/2;             cCloth(3)+lCloth/2];

% Define a rich trajectory for the upper corners
TrajU(:, 1:NT) = u_ini*ones(1,NT);

if (NExp == 8 || NExp == 9)
    TrajU(:,NT+1:3*NT) = TrajU(:,NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajU(:,3*NT+1:7*NT) = TrajU(:,3*NT)*ones(1,4*NT) + ...
                               [-linspace(0.0005,0.100,4*NT).*ones(2,1);
                                -linspace(0.0010,0.200,4*NT).*ones(2,1);
                                 zeros(2,4*NT)];
    TrajU(:,7*NT+1:9*NT) = TrajU(:,7*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajU(:,9*NT+1:11*NT) = TrajU(:,9*NT)*ones(1,2*NT) + ...
                                [linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 zeros(2,2*NT);];
    TrajU(:,11*NT+1:15*NT) = TrajU(:,11*NT)*ones(1,4*NT) + ...
                                [linspace(0.0010,0.200,4*NT).*ones(2,1);
                                 linspace(0.0010,0.200,4*NT).*ones(2,1);
                                -linspace(0.0005,0.100,4*NT).*ones(2,1)];
    TrajU(:,15*NT+1:17*NT) = TrajU(:,15*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0010,0.100,2*NT).*ones(2,1)];

    if NExp == 8
        % End directly (original traj)
        TrajU(:,17*NT+1:18*NT) = TrajU(:,17*NT)*ones(1,NT);
    elseif NExp == 9
        % Small rotation (Exp9 - YZ, same)
        TrajU(:,17*NT+1:19*NT) = TrajU(:,17*NT)*ones(1,2*NT) + ...
                                    [zeros(2,2*NT);
                                    -lCloth/2*sin(linspace(0.0005,pi/4,2*NT));
                                    -lCloth/2*sin(linspace(0.0005,pi/4,2*NT));
                                     lCloth/2*(1-cos(linspace(0.0005,pi/4,2*NT)));
                                     lCloth/2*(1-cos(linspace(0.0005,pi/4,2*NT)))];
        TrajU(:,19*NT+1:20*NT) = TrajU(:,19*NT)*ones(1,NT);
    end
elseif NExp == 10
    TrajU(:,NT+1:3*NT) = TrajU(:,NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajU(:,3*NT+1:5*NT) = TrajU(:,3*NT)*ones(1,2*NT) + ...
                               [zeros(2,2*NT);
                               -linspace(0.0010,0.200,2*NT).*ones(2,1);
                                zeros(2,2*NT)];
    TrajU(:,5*NT+1:6*NT) = TrajU(:,5*NT)*ones(1,NT);
    TrajU(:,6*NT+1:10*NT) = TrajU(:,6*NT)*ones(1,4*NT) + ...
                               [zeros(2,4*NT);
                                linspace(0.0010,0.400,4*NT).*ones(2,1);
                                zeros(2,4*NT)];
    TrajU(:,10*NT+1:11*NT) = TrajU(:,10*NT)*ones(1,NT);
    TrajU(:,11*NT+1:15*NT) = TrajU(:,11*NT)*ones(1,4*NT) + ...
                                    [zeros(2,4*NT);
                                    -0.2*(1-cos(linspace(0.0005,pi,4*NT)));
                                    -0.2*(1-cos(linspace(0.0005,pi,4*NT)));
                                    -0.2*sin(linspace(0.0005,pi,4*NT));
                                    -0.2*sin(linspace(0.0005,pi,4*NT))];
    TrajU(:,15*NT+1:16*NT) = TrajU(:,15*NT)*ones(1,NT);
    TrajU(:,16*NT+1:17*NT) = TrajU(:,16*NT)*ones(1,1*NT) + ...
                               [zeros(2,1*NT);
                                linspace(0.0020,0.200,1*NT).*ones(2,1);
                                linspace(0.0010,0.100,1*NT).*ones(2,1)];
elseif NExp == 11
    TrajU(:,NT+1:3*NT) = TrajU(:,NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajU(:,3*NT+1:5*NT) = TrajU(:,3*NT)*ones(1,2*NT) + ...
                                [linspace(0.0015,0.150,2*NT).*ones(2,1);
                                 zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajU(:,5*NT+1:9*NT) = TrajU(:,5*NT)*ones(1,4*NT) + ...
                               [-linspace(0.0020,0.200,4*NT).*ones(2,1);
                                -linspace(0.0020,0.200,4*NT).*ones(2,1);
                                 zeros(2,4*NT)];
    TrajU(:,9*NT+1:11*NT) = TrajU(:,9*NT)*ones(1,2*NT) + ...
                                   [-0.2*sin(linspace(0.0005,pi/2,2*NT));
                                    -0.2*sin(linspace(0.0005,pi/2,2*NT));
                                     0.2*(1-cos(linspace(0.0005,pi/2,2*NT)));
                                     0.2*(1-cos(linspace(0.0005,pi/2,2*NT)));
                                     zeros(2,2*NT);];
    TrajU(:,11*NT+1:13*NT) = TrajU(:,11*NT)*ones(1,2*NT) + ...
                                [linspace(0.0020,0.200,2*NT).*ones(2,1);
                                 zeros(2,2*NT);
                                 linspace(0.0010,0.100,2*NT).*ones(2,1)];
elseif NExp == 12
    TrajU(:,NT+1:3*NT) = TrajU(:,NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajU(:,3*NT+1:7*NT) = TrajU(:,3*NT)*ones(1,4*NT) + ...
                               [-linspace(0.0005,0.100,4*NT).*ones(2,1);
                                -linspace(0.0010,0.200,4*NT).*ones(2,1);
                                 zeros(2,4*NT)];
    TrajU(:,7*NT+1:9*NT) = TrajU(:,7*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajU(:,9*NT+1:11*NT) = TrajU(:,9*NT)*ones(1,2*NT) + ...
                                [linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 zeros(2,2*NT);];
    TrajU(:,11*NT+1:15*NT) = TrajU(:,11*NT)*ones(1,4*NT) + ...
                                [linspace(0.0010,0.200,4*NT).*ones(2,1);
                                 linspace(0.0010,0.200,4*NT).*ones(2,1);
                                -linspace(0.0005,0.100,4*NT).*ones(2,1)];
    TrajU(:,15*NT+1:17*NT) = TrajU(:,15*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0010,0.100,2*NT).*ones(2,1)];
    TrajU(:,17*NT+1:19*NT) = TrajU(:,17*NT)*ones(1,2*NT) + ...
                                [lCloth/2*(1-cos(linspace(0.0005,pi/4,2*NT)));
                                -lCloth/2*(1-cos(linspace(0.0005,pi/4,2*NT)));
                                 lCloth/2*sin(linspace(0.0005,pi/4,2*NT));
                                -lCloth/2*sin(linspace(0.0005,pi/4,2*NT));
                                 zeros(2,2*NT)];
    TrajU(:,19*NT+1:21*NT) = TrajU(:,19*NT:-1:17*NT+1);
    TrajU(:,21*NT+1:23*NT) = TrajU(:,21*NT)*ones(1,2*NT) + ...
                                [lCloth/2*(1-cos(linspace(0.0005,pi/4,2*NT)));
                                -lCloth/2*(1-cos(linspace(0.0005,pi/4,2*NT)));
                                -lCloth/2*sin(linspace(0.0005,pi/4,2*NT));
                                 lCloth/2*sin(linspace(0.0005,pi/4,2*NT));
                                 zeros(2,2*NT)];
elseif NExp == 13
    TrajU(:,NT+1:3*NT) = TrajU(:,NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajU(:,3*NT+1:7*NT) = TrajU(:,3*NT)*ones(1,4*NT) + ...
                               [-linspace(0.0005,0.100,4*NT).*ones(2,1);
                                -linspace(0.0010,0.200,4*NT).*ones(2,1);
                                 zeros(2,4*NT)];
    TrajU(:,7*NT+1:9*NT) = TrajU(:,7*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajU(:,9*NT+1:11*NT) = TrajU(:,9*NT)*ones(1,2*NT) + ...
                                [linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 zeros(2,2*NT);];
    TrajU(:,11*NT+1:15*NT) = TrajU(:,11*NT)*ones(1,4*NT) + ...
                                [linspace(0.0010,0.200,4*NT).*ones(2,1);
                                 linspace(0.0010,0.200,4*NT).*ones(2,1);
                                -linspace(0.0005,0.100,4*NT).*ones(2,1)];
    TrajU(:,15*NT+1:17*NT) = TrajU(:,15*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0010,0.100,2*NT).*ones(2,1)];
    TrajU(:,17*NT+1:19*NT) = TrajU(:,17*NT)*ones(1,2*NT) + ...
                                [lCloth/2*(1-cos(linspace(0.0005,pi/4,2*NT)));
                                -lCloth/2*(1-cos(linspace(0.0005,pi/4,2*NT)));
                                 zeros(2,2*NT)
                                 lCloth/2*sin(linspace(0.0005,pi/4,2*NT));
                                -lCloth/2*sin(linspace(0.0005,pi/4,2*NT))];
    TrajU(:,19*NT+1:21*NT) = TrajU(:,19*NT:-1:17*NT+1);
    TrajU(:,21*NT+1:23*NT) = TrajU(:,21*NT)*ones(1,2*NT) + ...
                                [lCloth/2*(1-cos(linspace(0.0005,pi/4,2*NT)));
                                -lCloth/2*(1-cos(linspace(0.0005,pi/4,2*NT)));
                                 zeros(2,2*NT);
                                -lCloth/2*sin(linspace(0.0005,pi/4,2*NT));
                                 lCloth/2*sin(linspace(0.0005,pi/4,2*NT))];
end

%% Get TrajTCP
cloth_x = TrajU([2 4 6],:) - TrajU([1 3 5],:);
cloth_y = [cloth_x(2,:);-cloth_x(1,:);zeros(1,size(cloth_x,2))];
cloth_z = cross(cloth_x,cloth_y);
cloth_x = permute(cloth_x./vecnorm(cloth_x), [1,3,2]);
cloth_y = permute(cloth_y./vecnorm(cloth_y), [1,3,2]);
cloth_z = permute(cloth_z./vecnorm(cloth_z), [1,3,2]);
Rcloth = [cloth_x cloth_y cloth_z];

Rtcp = [-cloth_y cloth_x cloth_z];
qtcp = rotm2quat(Rtcp); % Already in wxyz form

TrajTCP = zeros(7,size(TrajU,2));
TrajTCP(1,:) = mean(TrajU(1:2,:));
TrajTCP(2,:) = mean(TrajU(3:4,:));
TrajTCP(3,:) = mean(TrajU(5:6,:));
TrajTCP(4:7,:) = qtcp';

Ttcp = [Rtcp, permute(TrajTCP(1:3,:),[1,3,2]);
        repmat([0,0,0,1],1,1,size(Rtcp,3))];


%% Get initial joint positions
% WAM DH parameters
d1 = 0;
d3 = 0.55;
d5 = 0.3;
d7 = 0.09;

DH = [0 d1  0     -pi/2;
      0  0  0      pi/2;
      0 d3  0.045 -pi/2;
      0  0 -0.045  pi/2;
      0 d5  0     -pi/2;
      0  0  0      pi/2;
      0 d7  0        0];

% Create WAM Serial Link
L(1) = Link([DH(1,:), 0]);
L(2) = Link([DH(2,:), 0]);
L(3) = Link([DH(3,:), 0]);
L(4) = Link([DH(4,:), 0]);
L(5) = Link([DH(5,:), 0]);
L(6) = Link([DH(6,:), 0]);
L(7) = Link([DH(7,:), 0]);
wam = SerialLink(L, 'name', 'WAM [ALA]');

qref = [-pi/2 0.6 0 1.63, 0, 0.9, 0];
qini = wam.ikine(Ttcp(:,:,1), 'q0',qref);


%% Convert to traj for Cartesian Controller
TrajWAM = zeros(size(TrajTCP,2)+1,15);
TrajWAM(1,9:15) = TrajTCP(:,1)';
TrajWAM(2:end,9:15) = TrajTCP';
TrajWAM(:,1) = 0:Ts:Ts*size(TrajTCP,2);
TrajWAM(1,2:8) = qini;


%% Plot
fig1 = figure(1);
fig1.Color = [1,1,1];

plot3(TrajTCP(1,:), TrajTCP(2,:), TrajTCP(3,:))
axis equal
box on
grid on
hold on
plot3(TrajU(1,:), TrajU(3,:), TrajU(5,:))
plot3(TrajU(2,:), TrajU(4,:), TrajU(6,:))
hold off
set(gca, 'TickLabelInterpreter','latex');
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
zlabel('z','Interpreter','latex')

pov = [-50 30];
wamws = [-0.4 0.4 -1 0.2 -0.2 0.6];
wam.plot(qini, 'workspace', wamws, ...
               'notiles', 'noshadow', 'nobase', ...
               'jointdiam', 0.6, 'jointlen', 0.8, ...
               'lightpos', [0.4 0 1], 'fps', 30, ...
               'linkcolor', [1 0.6 0], 'view', pov, ...
               'jointcolor', [0.7 0 1], 'pjointcolor', [0.7 0 1]);


%% Animations
% Animation Frames
%{
anim = Animate(['TrajU_Mov',num2str(NExp)]);
qx = qini;
for tx=1:10:size(TrajU,2)
    Tx = [Rtcp(:,:,tx), TrajTCP(1:3,tx); 0 0 0 1];
    qx = wam.ikine(Tx, 'q0',qx);
    wam.plot(qx);
    anim.add();
end
%}
% Animation in plot
%{
qx = qini;
qxvect = zeros(size(TrajU,2),7);
for tx=1:size(TrajU,2)
    qx = wam.ikine(Ttcp(:,:,tx), 'q0',qx);
    qxvect(tx,:) = qx;
    if mod(tx,100)==0, fprintf(['IK timestep: ', num2str(tx), '\n']); end
end
wam.plot(qxvect, 'fps', 100);
%}


%% Save trajectories
if (savefiles)
    save(['TrajUs/TrajU_',num2str(NExp),'.mat'],'TrajU');
    save(['TrajTCPs/TrajTCP_',num2str(NExp),'.mat'],'TrajTCP');
    writematrix(TrajWAM,['TrajWAMs/TrajWAM_',num2str(NExp),'.csv']);
end




