function f = sourcefun()
% source term
Globals2D;

MN = 20;
[qx,qxw] = JacobiGQ(1-xorder, 0, MN); 
[qy,qyw] = JacobiGQ(1-yorder, 0, MN);
% 
tempx = zeros(Np,K); 
tempy = zeros(Np,K);
%%
% for k=1:K
%     for i=1:Np
%         etax = (x(i,k)+Ia)/2 + (x(i,k)-Ia)/2*qx;
%         etay = (y(i,k)+Ic)/2 + (y(i,k)-Ic)/2*qy;
%         tempx(i,k) = pi^2*((x(i,k) - Ia)/2)^(2-xorder)/gamma(2-xorder)*sin(pi*etax)'*qxw;
%         tempy(i,k) = pi^2*((y(i,k) - Ic)/2)^(2-yorder)/gamma(2-yorder)*sin(pi*etay)'*qyw;
%     end
% end
% f = tempx.*sin(pi*y) + sin(pi*x).*tempy - sin(pi*x).*sin(pi*y);
%%
for k=1:K
    for i=1:Np
        etax = (x(i,k)+Ia)/2 + (x(i,k)-Ia)/2*qx;
        etay = (y(i,k)+Ic)/2 + (y(i,k)-Ic)/2*qy;
%         tempx(i,k) = ((x(i,k)-Ia)/2).^(2-xorder)/gamma(2-xorder)*(2*pi*cos(pi*(etax+1).^2) - 4*pi^2*(etax+1).^2.*sin(pi*(etax+1).^2))'*qxw;
%         tempy(i,k) = ((y(i,k)-Ic)/2).^(2-yorder)/gamma(2-yorder)*(2*pi*cos(pi*(etay+1).^2) - 4*pi^2*(etay+1).^2.*sin(pi*(etay+1).^2))'*qyw;
        tempx(i,k) = ((x(i,k) - Ia)/2)^(2-xorder)/gamma(2-xorder)*(6*(etax.^2-1).^2 + 24*etax.^2.*(etax.^2-1))'*qxw;
        tempy(i,k) = ((y(i,k) - Ic)/2)^(2-yorder)/gamma(2-yorder)*(6*(etay.^2-1).^2 + 24*etay.^2.*(etay.^2-1))'*qyw;
    end
end
% f = - tempx.*sin(pi*(y+1).^2)  - sin(pi*(x+1).^2).*tempy - sin(pi*(x+1).^2).*sin(pi*(y+1).^2);
f = - tempx.*((y.^2-1).^3)  - ((x.^2-1).^3).*tempy - ((x.^2-1).^3).*((y.^2-1).^3);
%%
% f = zeros(Np,K); 
% for k=1:K
%     [xy] = transform(k); % Gauss integral points in the k-th triangle
%     w = trigauss(:,3);% Gauss weights in the standard triangle
%     for t = 1:length(xy)
%         tx = xy(t,1); ty = xy(t,2);
%         tmpx = pi^2/gamma(2-xorder)*((tx-Ia)/2)^(2-xorder)*sin(pi*((tx+Ia)/2+(tx-Ia)/2*qx))'*qxw;
%         tmpy = pi^2/gamma(2-yorder)*((ty-Ic)/2)^(2-yorder)*sin(pi*((ty+Ic)/2+(ty-Ic)/2*qy))'*qyw;
%         f(:,k) = f(:,k)  + ( tmpx*sin(pi*ty) + sin(pi*tx)*tmpy- sin(pi*tx)*sin(pi*ty))*LagrangeInt2DTRI(tx,ty,k)*w(t);
%     end
% end
% f = 4*invMassMatrix*f;
return;
