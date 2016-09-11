function f = fource()
% fource term of 2D advectio equation
Globals2D;

MN=20;
[qx,qxw] = JacobiGQ(-xorder, 0, MN);
[qy,qyw] = JacobiGQ(-yorder, 0, MN);

tempx = zeros(Np,K);
tempy = zeros(Np,K);

for k=1:K
    for i=1:Np
        etax = (x(i,k)+Ia)/2 + (x(i,k)-Ia)/2*qx;
        etay = (y(i,k)+Ic)/2 + (y(i,k)-Ic)/2*qy;
%         tempx(i,k) = pi^2*((x(i,k)-Ia)/2).^(2-xorder)/gamma(2-xorder)*cos(pi*etax)'*qxw;
        tempx(i,k) = ((x(i,k)-Ia)/2).^(1-xorder)/gamma(1-xorder)*cos(etax)'*qxw;
        tempy(i,k) = ((y(i,k)-Ic)/2).^(1-yorder)/gamma(1-yorder)*cos(etay)'*qyw;
    end
end
% f = sin(y).*tempx + sin(x).*cos(y) - sin(x).*sin(y);% x
f = sin(y).*tempx + sin(x).*tempy - sin(x).*sin(y);% both
% f = sin(y).*cos(x) + sin(x).*tempy - sin(x).*sin(y);% y
return;