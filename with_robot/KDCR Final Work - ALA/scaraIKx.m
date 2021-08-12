function qi = scaraIKx(xi, DH, flag, Eul)
%
% qi = scaraIKx(xi, DH, flag, Eul)
%
% Inverse Kinematics for a SCARA Manipulator:
%   Obtain joint variables vector qi
%   from a given end effector pose xi,
%   and a matrix of DH parameters.
%
% Use flag to choose one of the two possible solutions
%   it changes the sign of the sine of q2 (default at +1)
%
% As the orientation part on the pose vector xi can be represented
%   in different ways, input Eul has to be the Euler angle
%   where the orientation is saved (one variable angle needed for SCARA),
%   either 'phi' or 'psi'. If omitted, the function assumes that
%   xi(6) is the orientation of the end effector in the z axis, phiz
%   (same form as output xe from scaraFK).
%   With 'phi', phiz = x(4). With 'psi', phiz = -x(6).
%
% Example (computed xd and DH):
%   qd = scaraIKx(xd, DH, -1, 'psi')
%
% +-----------------------------------+
% | Adria Luque Acera (Alados5), 2020 |
% +-----------------------------------+

    if nargin < 4
        phiz = xi(6);
    elseif strcmp(Eul,'phi')
        phiz = xi(4);
    elseif strcmp(Eul, 'psi')
        phiz = -xi(6);
    else
        error(['Invalid Eul input.\n\nPlease specify either',...
               ' ''phi'' or ''psi'' on Euler angle input Eul.', ...
               '\nLeave input blank for default (z rotation)'],'');
    end

    if nargin < 3
        flag = 1;
    end
    
    pe = xi(1:3);
    
    th = DH(:,1);
    d = DH(:,2);
    a = DH(:,3);
    
    qi = zeros(1,4);
    
    qi(3) = d(1) + d(4) - d(3) - pe(3);
    
    c2 = (pe(1)^2+pe(2)^2 - a(1)^2 - a(2)^2)/(2*a(1)*a(2));
    s2 = sign(flag)*sqrt(1-c2^2);
    
    qi(2) = atan2(s2,c2);
    
    s1 = ((a(1)+a(2)*c2)*pe(2) - a(2)*s2*pe(1))/(pe(1)^2+pe(2)^2);
    c1 = ((a(1)+a(2)*c2)*pe(1) + a(2)*s2*pe(2))/(pe(1)^2+pe(2)^2);
    
    qi(1) = atan2(s1,c1);
    
    qi(4) = phiz - qi(1) - qi(2);
    
    qi = qi-th';
    qi([1:2,4]) = mod(qi([1:2,4])+pi,pi+pi)-pi;
    
end