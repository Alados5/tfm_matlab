clear; clc;

% Parameters
lCloth = 0.3;
cCloth = [0, -0.4, 0.1];
aCloth = 0;
NT = 50;
savefiles = 1;

% Ref Num: 
NRef = 16;

% Get initial positions 
r_ini = [cCloth(1)-lCloth/2*cos(aCloth); cCloth(1)+lCloth/2*cos(aCloth);
         cCloth(2)-lCloth/2*sin(aCloth); cCloth(2)+lCloth/2*sin(aCloth);
         cCloth(3)-lCloth/2;             cCloth(3)-lCloth/2];

% Define a trajectory for the lower corners
TrajR(:, 1:NT) = r_ini*ones(1,NT);

if (NRef == 10)
    TrajR(:,NT+1:3*NT) = TrajR(:,NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajR(:,3*NT+1:7*NT) = TrajR(:,3*NT)*ones(1,4*NT) + ...
                               [-linspace(0.0005,0.100,4*NT).*ones(2,1);
                                -linspace(0.0010,0.200,4*NT).*ones(2,1);
                                 zeros(2,4*NT)];
    TrajR(:,7*NT+1:9*NT) = TrajR(:,7*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajR(:,9*NT+1:11*NT) = TrajR(:,9*NT)*ones(1,2*NT) + ...
                                [linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 zeros(2,2*NT);];
    TrajR(:,11*NT+1:15*NT) = TrajR(:,11*NT)*ones(1,4*NT) + ...
                                [linspace(0.0010,0.200,4*NT).*ones(2,1);
                                 linspace(0.0010,0.200,4*NT).*ones(2,1);
                                -linspace(0.0005,0.100,4*NT).*ones(2,1)];
    TrajR(:,15*NT+1:17*NT) = TrajR(:,15*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0010,0.100,2*NT).*ones(2,1)];
    TrajR(:,17*NT+1:19*NT) = TrajR(:,17*NT)*ones(1,2*NT) + ...
                               [-linspace(0.0010,0.100,2*NT).*ones(2,1);
                                -linspace(0.0020,0.200,2*NT).*ones(2,1);
                                 zeros(2,2*NT)];
    TrajR(:,19*NT+1:20*NT) = TrajR(:,19*NT)*ones(1,1*NT) + ...
                               [-linspace(0.0010,0.100,1*NT).*ones(2,1);
                                 zeros(2,1*NT);
                                 zeros(2,1*NT)];
    TrajR(:,20*NT+1:22*NT) = TrajR(:,20*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                 linspace(0.0020,0.200,2*NT).*ones(2,1);
                                 linspace(0.0015,0.150,2*NT).*ones(2,1)];
    TrajR(:,22*NT+1:24*NT) = TrajR(:,22*NT)*ones(1,2*NT) + ...
                                [linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 zeros(2,2*NT);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajR(:,24*NT+1:25*NT) = TrajR(:,24*NT)*ones(1,1*NT);
        
elseif NRef == 11
    TrajR(:,NT+1:3*NT) = TrajR(:,NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajR(:,3*NT+1:7*NT) = TrajR(:,3*NT)*ones(1,4*NT) + ...
                               [-linspace(0.0005,0.100,4*NT).*ones(2,1);
                                -linspace(0.0010,0.200,4*NT).*ones(2,1);
                                 zeros(2,4*NT)];
    TrajR(:,7*NT+1:9*NT) = TrajR(:,7*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajR(:,9*NT+1:11*NT) = TrajR(:,9*NT)*ones(1,2*NT) + ...
                                [linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 zeros(2,2*NT);];
    TrajR(:,11*NT+1:15*NT) = TrajR(:,11*NT)*ones(1,4*NT) + ...
                                [linspace(0.0010,0.200,4*NT).*ones(2,1);
                                 linspace(0.0010,0.200,4*NT).*ones(2,1);
                                -linspace(0.0005,0.100,4*NT).*ones(2,1)];
    TrajR(:,15*NT+1:17*NT) = TrajR(:,15*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0010,0.100,2*NT).*ones(2,1)];
                            
    % Rotation in XY plane
    TrajR(:,17*NT+1:19*NT) = TrajR(:,17*NT)*ones(1,2*NT) + ...
                                [lCloth/2*(1-cos(linspace(0.0005,pi/6,2*NT)));
                                -lCloth/2*(1-cos(linspace(0.0005,pi/6,2*NT)));
                                 lCloth/2*sin(linspace(0.0005,pi/6,2*NT));
                                -lCloth/2*sin(linspace(0.0005,pi/6,2*NT));
                                 zeros(2,2*NT)];
    
    TrajR(:,19*NT+1:21*NT) = TrajR(:,19*NT)*ones(1,2*NT) + ...
                               [-linspace(0.0010,0.100,2*NT).*ones(2,1);
                                -linspace(0.0020,0.200,2*NT).*ones(2,1);
                                 zeros(2,2*NT)];
    TrajR(:,21*NT+1:22*NT) = TrajR(:,21*NT)*ones(1,1*NT) + ...
                               [-linspace(0.0010,0.100,1*NT).*ones(2,1);
                                 zeros(2,1*NT);
                                 zeros(2,1*NT)];
                             
    % Counter-Rotation in XY plane
    TrajR(:,22*NT+1:24*NT) = TrajR(:,19*NT:-1:17*NT+1) - TrajR(:,19*NT) + TrajR(:,22*NT);
                             
    TrajR(:,24*NT+1:26*NT) = TrajR(:,24*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                 linspace(0.0020,0.200,2*NT).*ones(2,1);
                                 linspace(0.0015,0.150,2*NT).*ones(2,1)];
    TrajR(:,26*NT+1:28*NT) = TrajR(:,26*NT)*ones(1,2*NT) + ...
                                [linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 zeros(2,2*NT);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajR(:,28*NT+1:29*NT) = TrajR(:,28*NT)*ones(1,1*NT);
    
elseif NRef == 12
    TrajR(:,NT+1:2*NT) = TrajR(:,1*NT)*ones(1,1*NT) + ...
                                [zeros(2,1*NT);
                                -linspace(0.0015,0.150,1*NT).*ones(2,1);
                                 zeros(2,1*NT)];
                             
    % Small circular motion (YZ, rise/lower wave form)
    TrajR(:,2*NT+1:3*NT) = TrajR(:,2*NT)*ones(1,1*NT) + ...
                                [zeros(2,1*NT);
                                -lCloth/2*sin(linspace(0.0010,pi/4,1*NT));
                                -lCloth/2*sin(linspace(0.0010,pi/4,1*NT));
                                 lCloth/2*(1-cos(linspace(0.0010,pi/4,1*NT)));
                                 lCloth/2*(1-cos(linspace(0.0010,pi/4,1*NT)))];
    TrajR(:,3*NT+1:4*NT) = TrajR(:,3*NT:-1:2*NT+1);
    TrajR(:,4*NT+1:5*NT) = TrajR(:,4*NT)*ones(1,1*NT) + ...
                                [zeros(2,1*NT);
                                 lCloth/2*sin(linspace(0.0010,pi/8,1*NT));
                                 lCloth/2*sin(linspace(0.0010,pi/8,1*NT));
                                 lCloth/2*(1-cos(linspace(0.0010,pi/8,1*NT)));
                                 lCloth/2*(1-cos(linspace(0.0010,pi/8,1*NT)))];
    TrajR(:,5*NT+1:6*NT) = TrajR(:,5*NT:-1:4*NT+1);
    
    TrajR(:,6*NT+1:8*NT) = TrajR(:,6*NT)*ones(1,2*NT) + ...
                               [-linspace(0.0010,0.100,2*NT).*ones(2,1);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0020,0.200,2*NT).*ones(2,1)];

	TrajR(:,8*NT+1:9*NT) = TrajR(:,8*NT)*ones(1,1*NT) + ...
                                [zeros(2,1*NT);
                                 linspace(0.0015,0.150,1*NT).*ones(2,1);
                                 zeros(2,1*NT)];
                             
    % Small circular motion (YZ, rise/lower wave form)
    TrajR(:,9*NT+1:10*NT) = TrajR(:,9*NT)*ones(1,1*NT) + ...
                                [zeros(2,1*NT);
                                -lCloth/2*sin(linspace(0.0010,-pi/4,1*NT));
                                -lCloth/2*sin(linspace(0.0010,-pi/4,1*NT));
                                 lCloth/2*(1-cos(linspace(0.0010,-pi/4,1*NT)));
                                 lCloth/2*(1-cos(linspace(0.0010,-pi/4,1*NT)))];
    TrajR(:,10*NT+1:11*NT) = TrajR(:,10*NT:-1:9*NT+1);
    TrajR(:,11*NT+1:12*NT) = TrajR(:,11*NT)*ones(1,1*NT) + ...
                                [zeros(2,1*NT);
                                 lCloth/2*sin(linspace(0.0010,-pi/8,1*NT));
                                 lCloth/2*sin(linspace(0.0010,-pi/8,1*NT));
                                 lCloth/2*(1-cos(linspace(0.0010,-pi/8,1*NT)));
                                 lCloth/2*(1-cos(linspace(0.0010,-pi/8,1*NT)))];
    TrajR(:,12*NT+1:13*NT) = TrajR(:,12*NT:-1:11*NT+1);
    
    TrajR(:,13*NT+1:14*NT) = TrajR(:,13*NT)*ones(1,1*NT) + ...
                                [zeros(2,1*NT);
                                 zeros(2,1*NT);
                                 linspace(0.0020,0.200,1*NT).*ones(2,1)];

    TrajR(:,14*NT+1:15*NT) = TrajR(:,14*NT)*ones(1,1*NT) + ...
                                [linspace(0.0010,0.100,1*NT).*ones(2,1);
                                 linspace(0.0005,0.050,1*NT).*ones(2,1);
                                 zeros(2,1*NT)];
    TrajR(:,15*NT+1:16*NT) = TrajR(:,15*NT)*ones(1,1*NT);
    
elseif NRef == 16
    TrajR(:,NT+1:3*NT) = TrajR(:,NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajR(:,3*NT+1:7*NT) = TrajR(:,3*NT)*ones(1,4*NT) + ...
                               [-linspace(0.0005,0.100,4*NT).*ones(2,1);
                                -linspace(0.0010,0.200,4*NT).*ones(2,1);
                                 zeros(2,4*NT)];
    TrajR(:,7*NT+1:9*NT) = TrajR(:,7*NT)*ones(1,2*NT) + ...
                                [zeros(2,2*NT);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1)];
    TrajR(:,9*NT+1:13*NT) = TrajR(:,9*NT)*ones(1,4*NT) + ...
                                [lCloth/2*(1-cos(linspace(0.0005,2*pi/3,4*NT)));
                                -lCloth/2*(1-cos(linspace(0.0005,2*pi/3,4*NT)));
                                -lCloth/2*sin(linspace(0.0005,2*pi/3,4*NT));
                                 lCloth/2*sin(linspace(0.0005,2*pi/3,4*NT));
                                 zeros(2,4*NT)];
    TrajR(:,13*NT+1:14*NT) = TrajR(:,13*NT)*ones(1,1*NT) + ...
                                [zeros(2,1*NT);
                                 zeros(2,1*NT);
                                -linspace(0.001,0.10,1*NT).*ones(2,1)];
    TrajR(:,14*NT+1:18*NT) = TrajR(:,14*NT)*ones(1,4*NT) + ...
                                [linspace(0.0005,0.10,4*NT).*ones(2,1);
                                 linspace(0.0005,0.10,4*NT).*ones(2,1);
                                 zeros(2,4*NT)];
    TrajR(:,18*NT+1:20*NT) = TrajR(:,18*NT)*ones(1,2*NT) + ...
                                [linspace(0.0005,0.050,2*NT).*ones(2,1);
                                -linspace(0.0005,0.050,2*NT).*ones(2,1);
                                 zeros(2,2*NT)];
                             
    TrajR(:,20*NT+1:21*NT) = TrajR(:,13*NT:-1:12*NT+1) - TrajR(:,13*NT) + TrajR(:,20*NT);
    TrajR(:,21*NT+1:22*NT) = TrajR(:,21*NT)*ones(1,1*NT) + ...
                                [zeros(2,1*NT);
                                 zeros(2,1*NT);
                                 linspace(0.0005,0.050,1*NT).*ones(2,1)];
    TrajR(:,22*NT+1:24*NT) = TrajR(:,22*NT)*ones(1,2*NT) + ...
                                [linspace(0.001,0.10,2*NT).*ones(2,1);
                                 zeros(2,2*NT);
                                 zeros(2,2*NT)];
    TrajR(:,24*NT+1:26*NT) = TrajR(:,24*NT)*ones(1,2*NT);
    
    
    
    
end


%% Reshape to ref file format
phiRef = TrajR';
phiL = phiRef(:,[1,3,5]);
phiR = phiRef(:,[2,4,6]);


%% Plot
fig1 = figure(1);
fig1.Color = [1,1,1];
plot3(phiL(:,1), phiL(:,2), phiL(:,3))
hold on
plot3(phiR(:,1), phiR(:,2), phiR(:,3))
hold off
axis equal
box on
grid on
set(gca, 'TickLabelInterpreter','latex');
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
zlabel('z','Interpreter','latex')

fig2 = figure(2);
fig2.Color = [1,1,1];
subplot(1,2,1)
plot(phiL)
grid on
set(gca, 'TickLabelInterpreter','latex');
subplot(1,2,2)
plot(phiR)
grid on
set(gca, 'TickLabelInterpreter','latex');




%% Save trajectories
if (savefiles==1)
    writematrix(phiL, ['trajectories/ref_',num2str(NRef),'L.csv'])
    writematrix(phiR, ['trajectories/ref_',num2str(NRef),'R.csv'])
end




