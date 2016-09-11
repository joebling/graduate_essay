function [p] = diffJacobiP(x,alpha,beta,N)
% first order derivative of normailzed Jacobi polynomial 
if (N==0)
    p=0;
else
    p=1/2*(N+alpha+beta+1)*JacobiP(x,alpha+1,beta+1,N-1);
end
return;
