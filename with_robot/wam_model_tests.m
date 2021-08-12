init_WAM;

q0 = [0 0 0 0 0 0 0];
q1 = [-0.00252 0.54735, -0.04479, 1.42524, -0.03874, 0.89746, -0.00462];
q2 = [0, 0.71, 0, 1.63, 0, 0.78, pi];
qH = [0.00256 -2.0115 0.00292 3.14977 -0.02127 -1.55826 -0.02723];

Rs = -eye(3); Rs(1)=1;
Ts = [Rs, [0.641514; 0; 0.15]; [0 0 0 1]];
qs =  wam.ikine(Ts, 'q0', q2);


%J  = wam.jacob0(q0);
%qtraj = jtraj([0 0 0 0], qr, linspace(0,1,60));

NTraj = 1;
TrajL = load(['../data/trajectories/ref_',num2str(NTraj),'L.csv']);
TrajR = load(['../data/trajectories/ref_',num2str(NTraj),'R.csv']);


TrajM_15 = load('./MPCtraj.csv');
TrajM = TrajM_15(:,9:11);

qt0 = qs;
Qt = zeros(size(TrajM,1),7);
for t=1:size(TrajM,1)
    pt = TrajM(t,:)';
    Tt = [Rs pt; 0 0 0 1];

    qt =  wam.ikine(Tt, 'q0',qt0, 'ilimit',500, 'rlimit',200, 'slimit',200, 'tol',1e-6);
    
    if isempty(qt)
        qt = qt0;
        disp('Warning!');
    end
    
    Qt(t,:) = qt;
    qt0 = qt;

end

%%

hf10 = figure(10);
hf10.Color = [1 1 1];

%set(hf10, 'menubar', 'none');

title('\textbf{WAM Robot Representation}', ...
      'Interpreter', 'latex','FontSize',14);
xlabel('X','Interpreter', 'latex','FontSize',12);
ylabel('Y','Interpreter', 'latex','FontSize',12);
zlabel('Z','Interpreter', 'latex','FontSize',12);
box on
lgr = legend('Location', 'SouthOutside', 'Orientation', 'horizontal', ...
             'Box', 'off', 'Interpreter', 'latex', 'FontSize', 14);

         
% Plot configuration
lt = 1;
pov = [30 36];
%{
wam.plot(qs, 'workspace', [-lt lt -lt lt -0.2*lt lt], ...
             'notiles', 'noshadow', 'nobase', ...
             'jointdiam', 0.6, 'jointlen', 0.8, ...
             'lightpos', [0.4 0 1], 'fps', 30, ...
             'linkcolor', [1 0.6 0], 'view', pov, ...
             'jointcolor', [0.7 0 1], 'pjointcolor', [0.7 0 1]);
%}
wam.plot(Qt, 'workspace', [-0.1*lt 0.8*lt -0.2*lt 0.5*lt -0.1*lt 0.7*lt], ...
             'notiles', 'noshadow', 'nobase', ...
             'jointdiam', 0.6, 'jointlen', 0.8, ...
             'lightpos', [0.4 0 1], 'fps', 30, ...
             'linkcolor', [1 0.6 0], 'view', pov, ...
             'jointcolor', [0.7 0 1], 'pjointcolor', [0.7 0 1]);
%
