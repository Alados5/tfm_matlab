function qi = scaraIK(Ti, DH, flag)
%
% qi = scaraIK(Ti, DH, flag)
%
% Inverse Kinematics for a SCARA Manipulator:
%   Obtain joint variables vector qi
%   from a given FK transform matrix Ti,
%   and a matrix of DH parameters
%
% Use flag to choose one of the two possible solutions
%   it changes the sign of the sine of q2 (default at +1)
%
% Example (computed Td and DH):
%   qd = scaraIK(Td, DH, -1)
%
% +-----------------------------------+
% | Adria Luque Acera (Alados5), 2020 |
% +-----------------------------------+

    if nargin < 3
        flag = 1;
    end
    
    if isa(Ti, 'double'), qi = zeros(1,4);
    else
        error(['Invalid data type for input: Ti',...
               '\nTi must be numerical'],'');
    end

    pe = Ti(1:3,4);
    ne = Ti(1:3,1);
    
    th = DH(:,1);
    d = DH(:,2);
    a = DH(:,3);
    
    qi(3) = d(1) + d(4) - d(3) - pe(3);
    
    c2 = (pe(1)^2+pe(2)^2 - a(1)^2 - a(2)^2)/(2*a(1)*a(2));
    s2 = sign(flag)*sqrt(1-c2^2);
    
    qi(2) = atan2(s2,c2);
    
    s1 = ((a(1)+a(2)*c2)*pe(2) - a(2)*s2*pe(1))/(pe(1)^2+pe(2)^2);
    c1 = ((a(1)+a(2)*c2)*pe(1) + a(2)*s2*pe(2))/(pe(1)^2+pe(2)^2);
    
    qi(1) = atan2(s1,c1);
    
    qi(4) = atan2(ne(2), ne(1)) - qi(1) - qi(2);
    
    qi = qi-th';
    qi([1:2,4]) = mod(qi([1:2,4])+pi,pi+pi)-pi;

end