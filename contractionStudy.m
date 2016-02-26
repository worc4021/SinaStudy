clear
close all
clc

Psi = [.2,5;0,.5];
opt = sdpsettings('solver','mosek-sdp');

Q = eye(2);
P = dlyap(Psi,Q);
rho = sdpvar(1);
Cons = Q>=rho*P;
optimize(Cons,-rho,opt);


Q1 = rand(2);
Q1 = Q1'*Q1;
P1 = dlyap(Psi,Q1);
rho1 = sdpvar(1);
Cons1 = Q1>=rho1*P1;
optimize(Cons1,-rho1,opt);

Q2 = Psi'*Psi;
P2 = dlyap(Psi,Q2);
rho2 = sdpvar(1);
Cons2 = Q2>=rho2*P2;
optimize(Cons2,-rho2,opt);


Q3 = diag([500,.5]);
P3 = dlyap(Psi,Q3);
rho3 = sdpvar(1);
Cons3 = Q3>=rho3*P3;
optimize(Cons3,-rho3,opt);


rho4 = 1 - (max(abs(eig(Psi))))^2;
P4 = sdpvar(2); 
alpha4 = sdpvar(1);
cons4 = [alpha4 >= 0,P4 >= alpha4*eye(2),P4 <= eye(2),(1-rho4)*P4-Psi'*P4*Psi >= 0];
info = optimize(cons4,-alpha4,opt);


[1-[value(rho),value(rho1),value(rho2),value(rho3),value(alpha4)]]