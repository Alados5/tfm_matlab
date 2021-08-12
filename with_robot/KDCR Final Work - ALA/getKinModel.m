function [T0n, q, An] = getKinModel(DH, sigma)
%
% [T0n, q, An] = getKinModel(DH, sigma)
%
% Get Kinematic Model (FK homogeneous transform) for a given
%   Denavit-Hartenberg parameter matrix DH and a joint type vector sigma.
%
% sigma can be either an array of 0 (revolute) / 1 (prismatic) or
%   a string with 'r' (revolute) / 'p' (prismatic)
%
% Matrix DH must have one row per joint, and the columns must be,
%   in order: theta, d, a, alpha. Parameters theta/d corresponding to
%   revolute/prismatic joints will be treated as offsets
%
% Example (after defining a 4-row matrix DH):
%   [T0n, q, An] = getKinModel(DH, 'rrpr')
%
% +-----------------------------------+
% | Adria Luque Acera (Alados5), 2020 |
% +-----------------------------------+

    % Number of Joints/Links
    n = length(sigma);
    if n ~= length(DH(:,1))
        error(['Incoherent number of joints.', ...
               '\nLength of sigma must be the number of rows of DH'],'');
    end

    % Joint types
    jtype = zeros(1,n);
    if ischar(sigma)
        sq = replace(sigma,'r','0');
        sq = replace(sq,'p','1');
        sq = num2cell(sq);
        
        for i=1:n
            jtype(i) = str2double(sq{i});
            if isnan(jtype(i))
                error(['Invalid sigma input.',...
                    '\n\nEnter ''r'' for revolute joints ', ...
                    'and ''p'' for prismatic joints', ...
                    '\nExample: ''rrpr'''],'');
            end
        end
        
    elseif isa(sigma, 'double')
        for i=1:n
            jtype(i) = logical(sigma(i));
        end
    else
        error(['Invalid data type for input: sigma', ...
               '\n\nType help getKinModel to see syntax'],'');
    end
    
    % Symbolic joint variables
    q = sym('q',[1 n]);
    assume(q,'real');

    % Symbolic DH Parameters
    syms thi di ai ali
    assume([thi di ai ali],'real')
    
    % General A matrix
    Ai = [cos(thi)  -sin(thi)*cos(ali)   sin(thi)*sin(ali)  ai*cos(thi);
          sin(thi)   cos(thi)*cos(ali)  -cos(thi)*sin(ali)  ai*sin(thi);
             0       sin(ali)            cos(ali)           di;
             0          0                   0                1];
    
    % DH Substitution and FK Transform
    T0n = eye(4);
    An = 0*sym('An',[4 4 n]);
    for i=1:n
        if jtype(i) == 0
            % Revolute joint: q+th d a alpha
            An(:,:,i) = subs(Ai, [thi, di,ai,ali], ...
                                 [q(i)+DH(i,1), DH(i,2:4)]);
        else
            % Prismatic joint: th q+d a alpha
            An(:,:,i) = subs(Ai, [thi, di, ai,ali], ...
                                 [DH(i,1), q(i)+DH(i,2), DH(i,3:4)]);
        end
        T0n = T0n*An(:,:,i);
    end
         
    T0n = simplify(T0n);
    An = simplify(An);
    
end