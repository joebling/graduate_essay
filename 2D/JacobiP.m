function [P]=JacobiP(x,alpha,beta,N)
% function [P]=JacobiP(x,alpha,beta,N)
% Purpose: Evaluate Jacobi Polynomial of type (alpha,beta)>-1 
%          at points x for order N and returns P[1:length(xp)]
%% Note:    They are normalized to be orthonormal

% Turn points into row if needed

xp=x; dims=size(xp);

if (dims(2)==1)
    xp=xp'; 
end

PL=zeros(2,length(xp));

% Initial values P_0(x) and P_1(x)
gamma0=2^(alpha+beta+1)/gamma(alpha+beta+2)*gamma(alpha+1)*gamma(beta+1);
PL(1,:)=1.0/sqrt(gamma0);
if (N==0) 
    P=(PL(1,:))';
     return;
end
gamma1=(alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
PL(2,:)=((alpha+beta+2)*xp/2+(alpha-beta)/2)/sqrt(gamma1);
if (N==1)
    P=(PL(2,:))'; 
     return; 
end

% Repeat value in recurrence
aold=2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

% Forward recurrence using the symmetry of the recurrence
for i=1:N-1
   h1=2*i+alpha+beta; 
   anew=2/(h1+2)*sqrt((i+1)*(i+alpha+beta+1)*(i+alpha+1)*(i+beta+1)/((h1+1)*(h1+3)));
   bnew=-(alpha^2-beta^2)/(h1*(h1+2));
   P=1/anew*(-aold*PL(1,:)+(xp-bnew).*PL(2,:));
   aold=anew;
   PL(1,:)=PL(2,:);
   PL(2,:)=P;
end
P=P';
return
