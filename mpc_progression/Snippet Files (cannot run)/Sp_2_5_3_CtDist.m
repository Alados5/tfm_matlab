% Constraint: Constant distance between upper corners
x_ctrl = x(cctrl,k+1);
g = [g; sum((x_ctrl([2,4,6]) - x_ctrl([1,3,5])).^2) - lCloth^2 ];
lbg = [lbg; -gbound];
ubg = [ubg;  gbound];