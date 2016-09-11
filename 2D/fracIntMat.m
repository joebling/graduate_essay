function [globalFDx, globalFDy] = fracIntMat( )
% global fractional integral matrix and 0<alpha<1

Globals2D;
MN = 10;
[qx,qxw] = JacobiGQ(-xorder,0,MN); 
[qy,qyw] = JacobiGQ(-yorder,0,MN);

globalFDx = zeros(K*Np,K*Np);
globalFDy = zeros(K*Np,K*Np);


for k=1:K % loop by elements 
    for i=1:Np % loop by interpolation points in every element
        tmpx = x(i,k);  tmpy = y(i,k);% loop by collocation points
%%
        [xdir,xreg] = findXdir(tmpx,tmpy);
        for j=1:length(xdir) % fractional matrix in x direction
            m = xdir(j); % the index of affected region in x direction
            vecx = zeros(1,Np); 
            if xreg(j,2)>100000 %%%%%
                a=xreg(j,1);
                if abs(tmpx-a)<NODETOL
                else
                    for kk=1:MN+1
                        vecx = vecx + ((tmpx-a)/2)^(1-xorder)*qxw(kk)*LagrangeInt2DTRI((tmpx+a)/2 + (tmpx-a)/2*qx(kk),tmpy,m)';% m is the index of element affected by collocation point
                    end
                end
            else %%%%%
                a=xreg(j,1);b=xreg(j,2);
                if abs(a-b)<NODETOL
                else
                    if abs(tmpx-b)<NODETOL
                        for kk=1:MN+1    
                            vecx = vecx + ((tmpx-a)/2)^(1-xorder)*qxw(kk)*LagrangeInt2DTRI((tmpx+a)/2 + (tmpx-a)/2*qx(kk),tmpy,m)';
                        end
                    else
                        for kk=1:MN+1
                            vecx = vecx + ((tmpx-a)/2)^(1-xorder)*qxw(kk)*LagrangeInt2DTRI((tmpx+a)/2 + (tmpx-a)/2*qx(kk),tmpy,m)' - ((tmpx-b)/2)^(1-xorder)*qxw(kk)*LagrangeInt2DTRI((tmpx+b)/2 + (tmpx-b)/2*qx(kk),tmpy,m)';
                        end
                    end
                end
            end
            globalFDx((k-1)*Np+i ,((m-1)*Np+1):m*Np) = 1/gamma(1-xorder)*vecx;
        end % end for m
 %%
        [ydir,yreg] = findYdir(tmpx,tmpy);
        for j=1:length(ydir)
            m=ydir(j); % the index of affected region in y direction
            vecy= zeros(1,Np); 
            if yreg(j,2)>100000 %%%%%
                c=yreg(j,1);
                if abs(tmpy - c)<NODETOL
                else
                    for kk=1:MN+1
                        vecy = vecy + ((tmpy-c)/2)^(1-yorder)*qyw(kk)*LagrangeInt2DTRI(tmpx,(tmpy+c)/2 + (tmpy-c)/2*qy(kk),m)';% m is index of element affected by collocation point
                    end
                end
            else %%%%%
                c=yreg(j,1);d=yreg(j,2);
                if abs(c-d)<NODETOL
                else
                    if abs(tmpy-d)<NODETOL
                        for kk=1:MN+1
                            vecy = vecy + ((tmpy-c)/2)^(1-yorder)*qyw(kk)*LagrangeInt2DTRI(tmpx,(tmpy+c)/2 + (tmpy-c)/2*qy(kk),m)';
                        end
                    else
                        for kk=1:MN+1
                            vecy = vecy + ((tmpy-c)/2)^(1-yorder)*qyw(kk)*LagrangeInt2DTRI(tmpx,(tmpy+c)/2 + (tmpy-c)/2*qy(kk),m)' - ((tmpy-d)/2)^(1-yorder)*qyw(kk)*LagrangeInt2DTRI(tmpx,(tmpy+d)/2 + (tmpy-d)/2*qy(kk),m)';
                        end
                    end
                end
            end % end if 
            globalFDy((k-1)*Np+i ,((m-1)*Np+1):(m*Np)) = 1/gamma(1-yorder)*vecy;
        end % end for
    end
end
return;