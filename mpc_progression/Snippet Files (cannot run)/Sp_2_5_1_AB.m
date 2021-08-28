A = [eye(3*N)         Ts*eye(3*N);
     (Ts/m)*Fp    eye(3*N)+(Ts/m)*Fv];
A(cctrl, 3*N+1:end) = 0;
A(cctrl+3*N, :) = 0;

B = zeros(2*3*N,6);
for i=1:length(cctrl)
    B(cctrl(i), i) = 1;
end