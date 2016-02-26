function dnu = normSim(t,nu,A,v1)

logA = logm(A);
At = expm(logA*t);
Am1t = expm(logA*(1-t));
M = Am1t'*Am1t;
n = v1'*M*v1;
dn = -v1'*(logA'*M+M*logA)*v1;
x = Am1t*v1./n;
dx = -logA*x-x*(dn/n^2);
N = At'*At;
dN = logA'*N+N*logA;
dnu = -diag(x)\diag(dx)*nu + dN*x + N*dx;