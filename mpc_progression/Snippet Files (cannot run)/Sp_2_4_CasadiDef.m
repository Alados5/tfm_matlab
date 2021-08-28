% Declare model variables
x = [SX.sym('pos',3*nr*nc,Hp+1);
     SX.sym('vel',3*nr*nc,Hp+1)];
u =  SX.sym('u',6,Hp);

% Initial parameters of the optimization problem
P  = SX.sym('P', 1+6, max(6*nc*nr, Hp+1)); 
x0 = P(1, :)';
Rp = P(1+(1:6), 1:Hp+1);

% Fill the state matrix so it depends only on x0 and u
x(:,1) = x0;
for k = 1:Hp
    x(:,k+1) = A*x(:,k) + B*u(:,k) + COM.dt*f_ext;
    
    % [Update objective function with x, u, Rp]
end
