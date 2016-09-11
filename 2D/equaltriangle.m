function [Nv, VX, VY, M, EToV] = equaltriangle(m)
% regular triangle mesh 
a = -1; b = 1;
c = -1; d = 1;
x = linspace(a,b,m+1);
VX = repmat(x,1,m+1);
VY = (reshape(VX,m+1,m+1))';
VY = VY(:)';
Nv = length(VX);

temp1 = zeros(m^2,3);
temp2 = temp1;
for i=1:m
    for j=1:m
        temp1((i-1)*m+j,:) = [(i-1)*(m+1)+j,(i-1)*(m+1)+j+1,i*(m+1)+j];
    end
end
for  i=1:m
    for j=1:m
        temp2((i-1)*m+j,:) = [(i)*(m+1)+j+1,(i)*(m+1)+j,(i-1)*(m+1)+j+1];
    end
end
EToV = zeros(2*m^2,3);
M = 2*m^2;
for k=1:m^2
    EToV(2*k-1,:) = temp1(k,:);
    EToV( 2*k ,:) = temp2(k,:);
end
M = 2*m^2;
return;

 