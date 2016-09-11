function [x,w]=JacobiGQ(alpha,beta,N)

% function [x,w]=JacobiGQ(alpha,beta,N)
% Purpose: Compute the N'th order Gauss quadrature points, x,
%          and weights, w, associated with the Jacobi Polynomial
%          of type (alpha,beta)> -1 (beta <> -0.5)

if (N==0)
    x(1)=-(alpha-beta)/(alpha+beta+2); 
    w(1)=2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+2); return;
end

% Form symmetric matrix from recurrence

%J=zeros(N+1);

h1=2*(1:N)+alpha+beta;
h2=2./h1.*sqrt((1:N).*((1:N)+alpha+beta).*((1:N)+alpha).*((1:N)+beta)./(h1-1)./(h1+1));
J=diag([-(alpha-beta)/(alpha+beta+2)  -(alpha^2-beta^2)./h1./(h1+2)])+...
    diag(h2,1)+diag(h2,-1);

if abs(alpha+beta+1)<eps
    J(1,2)=2*(1+alpha)*(1+beta)/((2+alpha+beta)^2*(alpha+beta+3));
    J(2,1)=J(1,2);
end
% Compute quadrature by eigenvalue solver 
[V,D]=eig(J); x=diag(D);
w=(V(1,:)').^2*2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+2);
return
