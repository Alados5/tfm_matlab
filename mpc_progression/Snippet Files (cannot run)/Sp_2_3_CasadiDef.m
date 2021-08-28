% Declare model variables
phi = SX.sym('phi',3*nr*nc);
dphi = SX.sym('dphi',3*nr*nc);
x = [phi; dphi];
u =  SX.sym('u',6);

% Mapping function f(x,u)->(x_next)
f = Function('f',{x,u},{A*x + B*u + COM.dt*f_ext}); 

% Parameters of the optimization problem: initial state and reference
P = SX.sym('P',6*nr*nc,Hp+1);

Xk = P(:,1);
for k = 1:Hp
    % Control variables
    Uk = SX.sym(['U_' num2str(k)],6);
    
    % Obtain the states for the next step
    Xk_next = f(Xk,Uk);
    Xk = Xk_next;
    
    % [Create symbolic artificial reference ra. Add it and Uk to Opt.v Vector]
    % [Update objective function with Xk_next, Uk, ra, P]
end