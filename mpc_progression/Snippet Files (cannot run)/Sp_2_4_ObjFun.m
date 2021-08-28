% Objective function
x_err = x(COM.coord_lc,k+1) - Rp(:,k+1);
objfun = objfun + x_err'*Q*x_err + u(:,k)'*R*u(:,k);