function [re] = transform(k)
% transform gauss integral pioints in standard triangle to specific
% triangle 
Globals2D;

B = [VX(EToV(k,2)) - VX(EToV(k,1)),VX(EToV(k,3)) - VX(EToV(k,1));
     VY(EToV(k,2)) - VY(EToV(k,1)),VY(EToV(k,3)) - VY(EToV(k,1))];
b = [VX(EToV(k,1));VY(EToV(k,1))];

re = (B*trigauss(:,1:2)' + b*ones(1,length(trigauss(:,1:2))))';
return;