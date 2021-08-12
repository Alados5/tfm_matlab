function qe = scaraIKiter(xd, Tq, JAq, Ts, q0, Alg, K, eStop, Eul, Timeout)
%
% qi = scaraIKiter(xd, Tq, JAq, Ts, q0, Alg, K, eStop, Eul, Timeout)
%
% Iterative Inverse Kinematics for a SCARA Manipulator:
%   From a given desired end effector pose xd,
%   the FK transform matrix Tq, the Analytic Jacobian JAq,
%   a sample time Ts, an initial joint vector q0,
%   the chosen algorithm Alg, a gain K, an error tolerance eStop,
%   and the corresponding Euler angle Eul,
%   get the joint configuration vector qe such that xd = FK(qe)
%
% A Timeout can be set indicating maximum iteration time
%   If omitted, the algorithm will stop after 60 seconds
%   if no solution is found, with an error message
%
% The general iterative method is valid for any robot, but the
%   FK calculation to get xe is done calling scaraFK, using also
%   the Eul parameter as one angle with the three possibilities
%   contemplated in that function.
%
% Pose vector xd, Analytic Jacobian JAq and Eul must be coherent,
%   or the function will not get correct results
%
% Example (computed Tq and JAq according to Eul):
%   [Td, xd] = scaraFK(Tq, [pi/2 -pi/3 0.3 -pi], 'phi');
%   qditer = scaraIKiter(xd, Tq, JAq, 1e-3, [0 0 0 0], 0, 100, 1e-4, 'phi');
%
% +-----------------------------------+
% | Adria Luque Acera (Alados5), 2020 |
% +-----------------------------------+

    if nargin < 10
        Timeout = 60;
    end
    if nargin < 9
        Eul = [];
    end
    if nargin < 8
        eStop = 1e-4;
    end
    if nargin < 7
        K = 100;
    end
    if nargin < 6
        Alg = 0;
    end
    if nargin < 5
        q0 = [0 0 0 0];
    end
    
    syms q1 q2 q3 q4
    assume([q1 q2 q3 q4], 'real');
    q = [q1 q2 q3 q4];
    
    qe = q0;
    [~, xe] = scaraFK(Tq,q0,Eul);
    e = xd-xe;
    eTh = max(abs(e));
    
    JAqinv = simplify(pinv(JAq));
    
    niter = 0;
    timer = tic;
    while eTh >= eStop

        if Alg == 0
            if qe(2) == 0 || qe(2) == pi
                JAinv = pinv(double(subs(JAq,q,qe)));
            else
                JAinv = double(subs(JAqinv, q, qe));
            end
            
            qdot = JAinv*K*e;
        else
            JA = double(subs(JAq, q, qe));
            qdot = JA'*K*e;
        end
    
        qe = qe + qdot'*Ts;
        qe([1:2,4]) = mod(qe([1:2,4])+pi,pi+pi)-pi;
        
        [~, xe] = scaraFK(Tq,qe,Eul);
        e = xd-xe;
        eTh = max(abs(e));
        niter = niter+1;
        
        disp(eTh);
        Titer = toc(timer);
        
        if Titer >= Timeout
            error(['Timeout error', '\n\nMaximum time reached', ...
                ' with no solution found.', ...
                '\nTry different input parameters.'],'');
        end
       
    end
    
end