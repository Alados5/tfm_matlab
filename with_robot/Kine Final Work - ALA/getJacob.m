function J = getJacob(An, sigma)
%
% J = getJacob(An, sigma)
%
% Get Geometric Jacobian for a given array of consecutive joint
%   transform matrices An (obtained using DH) and a joint type vector sigma.
%
% An must be a 3D array with 4x4 matrices in the first two dimensions,
%   corresponding to FK transforms between consecutive joints,
%   concatenated in the third dimension.
%   This matrix can contain symbolic variables corresponding
%   to joint variables q, then output J will be symbolic too.
%   The function works with An numerical too.
%
% sigma can be omitted only when An contains symbolic joint variables.
%   DO NOT omit it if it contains other symbolic variables or none.
%   If entered, it can be an array of 0 (revolute) / 1 (prismatic) or
%   a string with 'r' (revolute) / 'p' (prismatic)
%
% Example (in conjunction with getKinModel):
%   [T0n, q, An] = getKinModel(DH, 'rrpr');
%   J = getJacob(An, 'rrpr');
%
% +-----------------------------------+
% | Adria Luque Acera (Alados5), 2020 |
% +-----------------------------------+

    n = length(An(1,1,:));
    jtype = zeros(1,n);
    
    if nargin < 2
        if isa(An, 'double')
            error(['Joint types unknown for numerical transforms', ...
                   '\n\nsigma can only be omitted with symbolic An'], '');
        elseif isa(An, 'sym')
            
            % Joint types deduction
            for i=1:n
                jtype(i) = ~isempty(symvar(An(3,4,i)));
            end
            
        else
            error(['Invalid data type for input: An', ...
                   '\n\nAn must be numerical or symbolic An(q)'], '');
        end
    else
        if n ~= length(sigma)
        error(['Incoherent number of joints.', ...
               '\nLength of sigma and 3rd dimension of An must be the same'],'');
        end
        
        % Joint types from input
        if ischar(sigma)
            sq = replace(sigma,'r','0');
            sq = replace(sq,'p','1');
            sq = num2cell(sq);
            
            for i=1:n
                jtype(i) = str2double(sq{i});
                if isnan(jtype(i)), error('Invalid sigma input'); end
            end
            
        elseif isa(sigma, 'double')
            for i=1:n
                jtype(i) = logical(sigma(i));
            end
        else
            error('Invalid data type for input: sigma');
        end
        
    end
    
    p = 0*sym('p',[4,1,n+1]);
    z = 0*sym('z',[3,1,n+1]);
    
    p0 = [0 0 0 1]';
    z0 = [0 0 1]';
    
    % Index i corresponds to p and z of i-1
    p(:,:,1) = p0;
    z(:,:,1) = z0;
    Tp = eye(4);
    for i=1:n
        Ai = An(:,:,i);
        Tp = Tp*Ai;
        p(:,:,i+1) = Tp*p0;
        z(:,:,i+1) = Tp(1:3,1:3)*z0;
    end

    J = 0*sym('J',[6,n]);
    for i=1:n
        if jtype(i) == 0
            % Revolute joint
            J(1:3,i) = cross(z(:,:,i), p(1:3,:,n+1)-p(1:3,:,i));
            J(4:6,i) = z(:,:,i);
        else
            % Prismatic joint
            J(1:3,i) = z(:,:,i);
            J(4:6,i) = [0,0,0]';
        end
    end
    
    J = simplify(J);
    
    if isa(An,'double')
        J = double(J);
    end

end