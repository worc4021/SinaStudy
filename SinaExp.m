clear
close all
clc


h = 1e-1;
dissipation = 1e-2;
tE = 10;
x0 = [1;1];

A = [0,1;-1,0];
C = zeros(2);
C(1,1) = -dissipation;

PSI = A+C;

f = @(t,x)PSI*x;

[t,Y] = ode45(f,linspace(0,tE,100),x0);

Ae = eye(2)+h*PSI;
Ai = (eye(2)-h*PSI)\eye(2);
As = ([0,1;1+h*dissipation,-h]\[-h,1;1,0]);
Aexp = expm(PSI*h);

N = tE/h;
Xe = [x0,zeros(2,N)];
Xi = [x0,zeros(2,N)];
Xs = [x0,zeros(2,N)];
XExp = [x0,zeros(2,N)];
Tt = linspace(0,tE,N+1);


for i = 1:N
    Xe(:,i+1) = Ae*Xe(:,i);
    Xi(:,i+1) = Ai*Xi(:,i);
    Xs(:,i+1) = As*Xs(:,i);
    XExp(:,i+1) = Aexp*XExp(:,i);
end
figure(1)
plot(t,Y(:,1))
hold on
plot(Tt,Xe(1,:))
plot(Tt,Xi(1,:))
plot(Tt,Xs(1,:))
plot(Tt,XExp(1,:))
legend('Runge-Kutta','Explicit-Euler','Implicit-Euler','Simplectic',...
    'Exact Discretisation');
hold off

F = [1,0;-1,0]*1/10;
f = ones(2,1);
LambdaE = F;
lambdaE = f;
LambdaI = F;
lambdaI = f;
LambdaS = F;
lambdaS = f;
LambdaExp = F;
lambdaExp = f;


LambdaENext = [];
lambdaENext = [];

iter = 1;
iterMax = 50;
while and(~isContained(LambdaE,lambdaE,LambdaENext,lambdaENext),iter<iterMax)
    
    
    if iter~=1
        LambdaE = LambdaENext;
        lambdaE = lambdaENext;
    end
    
    LambdaENext=LambdaE*Ae;
    lambdaENext=ones(size(lambdaE));
    
    [LambdaENext,lambdaENext] = inequalityReduction([LambdaE;LambdaENext],...
        [lambdaE;lambdaENext]);
    iter = iter+1;

end

eLambda = LambdaENext;
elambda = lambdaENext;

LambdaINext = [];
lambdaINext = [];

iter = 1;
while and(~isContained(LambdaI,lambdaI,LambdaINext,lambdaINext),iter<iterMax)
    
    
    if iter~=1
        LambdaI = LambdaINext;
        lambdaI = lambdaINext;
    end
    
    LambdaINext=LambdaI*Ai;
    lambdaINext=ones(size(lambdaI));
    
    [LambdaINext,lambdaINext] = inequalityReduction([LambdaI;LambdaINext],...
        [lambdaI;lambdaINext]);
    iter = iter+1;

end

iLambda = LambdaINext;
ilambda = lambdaINext;


LambdaSNext = [];
lambdaSNext = [];

iter = 1;
iterMax = 50;
while and(~isContained(LambdaS,lambdaS,LambdaSNext,lambdaSNext),iter<iterMax)
    
    
    if iter~=1
        LambdaS = LambdaSNext;
        lambdaS = lambdaSNext;
    end
    
    LambdaSNext=LambdaS*As;
    lambdaSNext=ones(size(lambdaS));
    
    [LambdaSNext,lambdaSNext] = inequalityReduction([LambdaS;LambdaSNext],...
        [lambdaS;lambdaSNext]);
    iter = iter+1;

end

sLambda = LambdaSNext;
slambda = lambdaSNext;

LambdaExpNext = [];
lambdaExpNext = [];

iter = 1;
while and(~isContained(LambdaExp,lambdaExp,LambdaExpNext,lambdaExpNext),iter<iterMax)
    
    
    if iter~=1
        LambdaExp = LambdaExpNext;
        lambdaExp = lambdaExpNext;
    end
    
    LambdaExpNext=LambdaExp*Aexp;
    lambdaExpNext=ones(size(lambdaExp));
    
    [LambdaExpNext,lambdaExpNext] = inequalityReduction([LambdaExp;LambdaExpNext],...
        [lambdaExp;lambdaExpNext]);
    iter = iter+1;

end

expLambda = LambdaExpNext;
explambda = lambdaExpNext;



P = sdpvar(2);
CONS = [P>=1e-5*eye(2),PSI'*P+P*PSI<=-1e-5*eye(2)];
optimize(CONS,trace(P),sdpsettings('solver','mosek-sdp','verbose',0));

[V,~] = eig(value(P));

figure(2)
plot(Polyhedron(eLambda,elambda),'alpha',.3,'color','r',...
    Polyhedron(iLambda,ilambda),'alpha',.3,'color','b',...
    Polyhedron(sLambda,slambda),'alpha',.3,'color','g',...
    Polyhedron(expLambda,explambda),'alpha',.3,'color','y')
hold on

Z = V'*[10*sin(linspace(0,2*pi,250));...
     10*cos(linspace(0,2*pi,250))];

plot(Z(1,:),Z(2,:),'LineWidth',3)
hold off
