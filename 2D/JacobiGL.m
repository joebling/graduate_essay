function [x,w]=JacobiGL(alpha,beta,N)

% function [x]=JacobiGL(alpha,beta,N)
% Purpose: Compute the N'th order Gauss Lobatto quadrature
%          points, x, associated with the Jacobi polynomial, 
%          of type (alpha,beta) > -1 (beta <> 0.5).

x=zeros(N+1,1);
w=zeros(N+1,1);
w1=2^(alpha+beta+1)*gamma(beta+1)*gamma(beta+2)*gamma(N)*gamma(N+alpha+1)/gamma(N+beta+1)/gamma(N+alpha+beta+2);
wnplus1=2^(alpha+beta+1)*gamma(alpha+1)*gamma(alpha+2)*gamma(N)*gamma(N+beta+1)/gamma(N+alpha+1)/gamma(N+alpha+beta+2);

if (N==1)
    x(1)=-1.0; x(2)=1.0; 
    w(1)=w1; w(2)=wnplus1;   return; 
end

[xint,]=JacobiGQ(alpha+1,beta+1,N-2);
x=[-1,xint',1]';

[P]=JacobiP(xint,alpha,beta,N);
w(1)=w1;
w(N+1)=wnplus1;
w(2:N)=(2*N+alpha+beta+1)/(N*(N+alpha+beta+1))*(1-xint.^2)./P.^2./((4*(N+alpha)*(N+beta)+(alpha-beta)^2)/(2*N+alpha+beta)^2-xint.^2);
return
