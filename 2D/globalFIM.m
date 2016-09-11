function [globalFIx,globalFIy]=globalFIM()
% form the fractional stiffness matrix in a direct way
Globals2D;

globalFIx = zeros(K*Np,K*Np);
globalFIy = zeros(K*Np,K*Np);

MN = 10;
[qx,qxw] = JacobiGQ(1-xorder,0,MN);
[qy,qyw] = JacobiGQ(1-yorder,0,MN);

for k=1:K
    [xy] = transform(k); % Gauss integral points in the k-th triangle
%     det(B)
    w = trigauss(:,3);% Gauss weights in the standard triangle
    for t = 1:length(xy)
        tmpx = xy(t,1); tmpy = xy(t,2);% Gauss integral points in the specific triangle
        [xdir,xreg] = findXdir(tmpx,tmpy);% find the index and region influenced by Gauss integral points of the k-th element
       %% form fractional integral mass matrix in x direction
        for j=1:length(xdir)
            m = xdir(j);% the mark of affected triangle 
            matx = zeros(Np,Np);%mass matrix formed by affected domain
            if xreg(j,2)>100000 % the right end of the interval is Inf
                a=xreg(j,1); % the horizontal ordinate of the initial point 
                if abs(tmpx-a)<NODETOL %
                else
                    vecx = zeros(1,Np);
                    for kk =1:length(qx)
                        vecx = vecx + 1/gamma(2-xorder)*((tmpx-a)/2)^(2-xorder)*qxw(kk)*LagrangeInt2DTRI((tmpx+a)/2 + (tmpx-a)/2*qx(kk),tmpy,m)';
                    end
                    matx = matx +  LagrangeInt2DTRI(tmpx,tmpy,k)*vecx*w(t);
                end % end if 2
            else %%%%%
                a=xreg(j,1);b=xreg(j,2);
                if abs(a-b)<NODETOL %
                else
                    if abs(tmpx-b)<NODETOL% 
                        vecx = zeros(1,Np);
                        for kk =1:length(qx)
                            vecx = vecx + 1/gamma(2-xorder)*((tmpx-a)/2)^(2-xorder)*qxw(kk)*LagrangeInt2DTRI((tmpx+a)/2 + (tmpx-a)/2*qx(kk),tmpy,m)';
                        end
                        matx = matx +  LagrangeInt2DTRI(tmpx,tmpy,k)*vecx*w(t);
                    else
                        vecx = zeros(1,Np);
                        for kk =1:length(qx)
                            vecx = vecx + 1/gamma(2-xorder)*(((tmpx-a)/2)^(2-xorder)*qxw(kk)*LagrangeInt2DTRI((tmpx+a)/2 + (tmpx-a)/2*qx(kk),tmpy,m)'-...
                                ((tmpx-b)/2)^(2-xorder)*qxw(kk)*LagrangeInt2DTRI((tmpx+b)/2 + (tmpx-b)/2*qx(kk),tmpy,m)');
                        end
                        matx = matx +  LagrangeInt2DTRI(tmpx,tmpy,k)*vecx*w(t);
                    end
                end % end if 3
            end %end if 1 
            globalFIx((k-1)*Np+1:k*Np ,(m-1)*Np+1:m*Np) = globalFIx((k-1)*Np+1:k*Np ,(m-1)*Np+1:m*Np) + 4*matx;
        end % end for
       %% end for fractional integral matrix in x direction
       %% form fractional integral mass matrix in y direction
        [ydir,yreg] = findYdir(tmpx,tmpy);
       for j=1:length(ydir)
           m = ydir(j);
           maty = zeros(Np,Np);
           if yreg(j,2)>100000 %%%%%
                c=yreg(j,1);
                if abs(tmpy - c)<NODETOL
                else
                    vecy = zeros(1,Np);
                    for kk=1:length(qy)
                        vecy = vecy + 1/gamma(2-yorder)*((tmpy-c)/2)^(2-yorder)*qyw(kk)*LagrangeInt2DTRI(tmpx,(tmpy+c)/2 + (tmpy-c)/2*qy(kk),m)';
                    end
                    maty = maty +  LagrangeInt2DTRI(tmpx,tmpy,k)*vecy*w(t);
                end
            else %%%%%
                c=yreg(j,1);d=yreg(j,2);
                if abs(c-d)<NODETOL
                else
                    if abs(tmpy-d)<NODETOL
                        vecy = zeros(1,Np);
                        for kk=1:length(qy)
                            vecy = vecy + 1/gamma(2-yorder)*((tmpy-c)/2)^(2-yorder)*qyw(kk)*LagrangeInt2DTRI(tmpx,(tmpy+c)/2 + (tmpy-c)/2*qy(kk),m)';
                        end
                        maty = maty +  LagrangeInt2DTRI(tmpx,tmpy,k)*vecy*w(t);
                    else
                        vecy = zeros(1,Np);
                        for kk=1:length(qy)
                            vecy = vecy + 1/gamma(2-yorder)*(((tmpy-c)/2)^(2-yorder)*qyw(kk)*LagrangeInt2DTRI(tmpx,(tmpy+c)/2 + (tmpy-c)/2*qy(kk),m)' - ...
                                ((tmpy-d)/2)^(2-yorder)*qyw(kk)*LagrangeInt2DTRI(tmpx,(tmpy+d)/2 + (tmpy-d)/2*qy(kk),m)');
                        end
                        maty = maty +  LagrangeInt2DTRI(tmpx,tmpy,k)*vecy*w(t);
                    end
                end
            end % end if
            globalFIy((k-1)*Np+1:k*Np ,(m-1)*Np+1:m*Np) = globalFIy((k-1)*Np+1:k*Np ,(m-1)*Np+1:m*Np) + 4*maty;
       end
    end % end t 
end
return;