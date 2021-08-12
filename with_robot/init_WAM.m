% WAM DH parameters
d1 = 0;
d3 = 0.56;
d5 = 0.31;
d7 = 0.05;

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