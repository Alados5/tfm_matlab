function [Te, xe] = scaraFK(Tq, qi, Eul)
%
% [Te, xe] = scaraFK(Tq, qi)
%
% Get Forward Kinematics Transform Matrix Te and end effector
%   pose vector xe = [px py pz phix phiy phiz]' for a SCARA Manipulator
%   from a given symbolic Transform matrix Tq and configuration qi.
%   If Eul is used, then xe = [px py pz phi theta psi]', so that
%   always theta = pi, and the orientation phiz 
%   is saved in the input angle Eul ('phi' or 'psi').
%   If Eul is omitted, then xe = [px py pz 0 0 phiz]' by default.
%
% The function substitutes joint variables in Tq by values in qi
%   and extracts coordinates x y z from the fourth column
%   (valid for any robot)
%   But also gets the phiz orientation using formulas specific
%   for the SCARA Manipulator, and assigns it to the specified
%   angle with the correct formulation.
%
% Examples (with computed Tq):
%   [Te, xe] = scaraFK(Tq, [pi/2 -pi/2 0.2 pi/2])
%   [Te, xe] = scaraFK(Tq, [pi/2 -pi/2 0.2 pi/2], 'psi')
%
% +-----------------------------------+
% | Adria Luque Acera (Alados5), 2020 |
% +-----------------------------------+

    if nargin < 3 || isempty(Eul)
        Eul = 0;
    elseif strcmp(Eul,'phi')
        Eul = 1;
    elseif strcmp(Eul, 'psi')
        Eul = 3;
    else
        error(['Invalid Eul input.\n\nPlease specify either',...
               ' ''phi'' or ''psi'' on Euler angle input Eul.', ...
               '\nLeave input blank for default (z rotation)'],'');
    end

    q = symvar(Tq);
    assume(q,'real')
    
    Te = double(subs(Tq, q, qi));
    pe = Te(1:3,4);
    
    % SCARA-Specific, geometric rotation:
    phiz = atan2(Te(2,1), Te(1,1));
    
    % Assigned to requested angle:
    if Eul == 0
        phie = [0;0;phiz];
    elseif Eul == 1
        phie = [phiz;pi;0];
    else
        phie = [0;pi;-phiz];
    end

    % Final pose vector:
    xe = [pe; phie];
    
end