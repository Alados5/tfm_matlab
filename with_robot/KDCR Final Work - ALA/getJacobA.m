function JA = getJacobA(Tq, Eul)
%
% J = getJacobA(Tq, Eul)
%
% Get Analytic Jacobian for a given FK transform matrix Tq
%   depending on symbolic joint variables.
%
% Eul describes the chosen Euler angles: [phi;theta;psi]
%   If omitted, the function calculates with the default
%   ZYZ expressions, which might not work.
%
% Example (in conjunction with getKinModel) with SCARA DH:
%   [T04, q, Ai4] = getKinModel(DH, 'rrpr');
%   JA = getJacobA(T04, [0;0;q(1)+q(2)+q(4)]);
%
% +-----------------------------------+
% | Adria Luque Acera (Alados5), 2020 |
% +-----------------------------------+

    TPhie = Tq(1:3, 1:3);

    if nargin < 2
        phi = atan2(TPhie(2,3), TPhie(1,3));
        theta = atan2(sqrt(TPhie(1,3)^2 + TPhie(2,3)^2), TPhie(3,3));
        psi = atan2(TPhie(3,2), -TPhie(3,1));
    else
        phi = Eul(1);
        theta = Eul(2);
        psi = Eul(3);
    end

    Pe = Tq(1:3,4);
    Phie = [phi;theta;psi];
    
    q = symvar(Tq);
    n = length(q);
    
    JA = 0*sym('JA',[6 n]);
    for i=1:n
        JA(1:3,i) = diff(Pe,   q(i));
        JA(4:6,i) = diff(Phie, q(i));
    end
    JA = simplify(JA);
    
end




