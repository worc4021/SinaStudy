clear
close all
clc

tol = 1e-16;
h = 1e-1;
dissipation = 1e-1;
tE = 10;
x0 = [1;1];

wMax = dissipation/50;

A = [0,1;-1,0];
C = zeros(2);
C(1,1) = -dissipation;

PSI = A+C;
d = [0;1];

Ae = eye(2)+h*PSI;
de = h*d;

Ai = (eye(2)-h*PSI)\eye(2);
di = h*(eye(2)-h*PSI)\d;

As = [0,1;1+h*dissipation,-h]\[-h,1;1,0];
ds = h*[0,1;1+h*dissipation,-h]\d;

Aexp = expm(PSI*h);
dexp = PSI\(Aexp-eye(2))*d;

Lambda = [1,0;-1,0]*1/10;
lambda = ones(2,1);

kMax = 5000;

CompI = zeros(2,kMax);
CompS = zeros(2,kMax);
CompX = zeros(2,kMax);

CompI(:,1) = abs(Lambda*di);
CompS(:,1) = abs(Lambda*ds);
CompX(:,1) = abs(Lambda*dexp);

k = 1;

AkI = Ai*di;
AkS = As*ds;
AkX = Aexp*dexp;

while k<kMax && (max(CompI(:,k))>tol || max(CompS(:,k))>tol || max(CompX(:,k))>tol)
    
    CompI(:,k+1) = abs(Lambda*AkI);
    CompS(:,k+1) = abs(Lambda*AkS);
    CompX(:,k+1) = abs(Lambda*AkX);
    
    AkI = Ai*AkI;
    AkS = As*AkS;
    AkX = Aexp*AkX;
    k = k+1;
end

CandI = max(sum(CompI,2));
CandS = max(sum(CompS,2));
CandX = max(sum(CompX,2));


LambdaI = [Lambda,zeros(size(lambda));
            zeros(2),[CandI;-1]];
lambdaI = [lambda;1;0];
LambdaINext = [];
lambdaINext = [];

iterI = 1;
iterIMax = 250;
while and(~isContained(LambdaI,lambdaI,LambdaINext,lambdaINext),iterI<iterIMax)
    
    
    if iterI~=1
        LambdaI = LambdaINext;
        lambdaI = lambdaINext;
    end
    
    LambdaINext=[LambdaI(:,1:2)*Ai,abs(LambdaI(:,1:2)*di)+LambdaI(:,3)];
    lambdaINext=lambdaI;
    
    [LambdaINext,lambdaINext] = bigReduce([LambdaI;LambdaINext],...
        [lambdaI;lambdaINext]);
    iterI = iterI+1;

end

iLambda = LambdaINext;
ilambda = lambdaINext;

LambdaS = [Lambda,zeros(size(lambda));
            zeros(2),[CandS;-1]];
lambdaS = [lambda;1;0];
LambdaSNext = [];
lambdaSNext = [];

iterS = 1;
iterSMax = 500;
while and(~isContained(LambdaS,lambdaS,LambdaSNext,lambdaSNext),iterS<iterSMax)
    
    
    if iterS~=1
        LambdaS = LambdaSNext;
        lambdaS = lambdaSNext;
    end
    
    LambdaSNext=[LambdaS(:,1:2)*As,abs(LambdaS(:,1:2)*ds)+LambdaS(:,3)];
    lambdaSNext=lambdaS;
    
    [LambdaSNext,lambdaSNext] = bigReduce([LambdaS;LambdaSNext],...
        [lambdaS;lambdaSNext]);
    iterS = iterS+1;

end

sLambda = LambdaSNext;
slambda = lambdaSNext;


LambdaEX = [Lambda,zeros(size(lambda));
            zeros(2),[CandX;-1]];
lambdaEX = [lambda;1;0];
LambdaEXNext = [];
lambdaEXNext = [];

iterEX = 1;
iterEXMax = 500;
while and(~isContained(LambdaEX,lambdaEX,LambdaEXNext,lambdaEXNext),iterEX<iterEXMax)
    
    
    if iterEX~=1
        LambdaEX = LambdaEXNext;
        lambdaEX = lambdaEXNext;
    end
    
    LambdaEXNext=[LambdaEX(:,1:2)*Aexp,abs(LambdaEX(:,1:2)*dexp)+LambdaEX(:,3)];
    lambdaEXNext=lambdaEX;
    
    [LambdaEXNext,lambdaEXNext] = bigReduce([LambdaEX;LambdaEXNext],...
        [lambdaEX;lambdaEXNext]);
    iterEX = iterEX+1;

end

exLambda = LambdaEXNext;
exlambda = lambdaEXNext;