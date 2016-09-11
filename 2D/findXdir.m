function [xdir,xreg]=findXdir(px,py)
Globals2D;
xdir=zeros(K,1);
xreg=Inf(K,2);
i=1;
for k=1:K
    miny = min( VY(EToV(k,:))); maxy = max( VY(EToV(k,:)));
    if  py>=miny && py<=maxy
        re= [];
        vx = VX(EToV(k,:));  vy = VY(EToV(k,:));
        if (vy(1) - py)*(vy(2) - py)<=0
%             if abs(vy(2)-vy(1))<=NODETOL
%                 re = [re,px];
%             else
                re = [re,(vx(2) - vx(1))/(vy(2)-vy(1))*(py-vy(1)) + vx(1)];
%             end
        end
        if (vy(1) - py)*(vy(3) - py)<=0
%             if abs(vy(3)-vy(1))<=NODETOL
%                 re = [re,px];
%             else
                re = [re,(vx(3) - vx(1))/(vy(3)-vy(1))*(py-vy(1)) + vx(1)];
%             end
        end
        if (vy(3) - py)*(vy(2) - py)<=0
%             if abs(vy(3)-vy(2))<=NODETOL
%                 re = [re,px];
%             else
                re = [re,(vx(3) - vx(2))/(vy(3)-vy(2))*(py-vy(2)) + vx(2)];
%             end
        end
        re = re(re<=px+NODETOL);
       %%
        if length(re)==1
            xdir(i)=k;
            xreg(i,1) = re; i=i+1; 
        elseif length(re)==2
            re = sort(re);  xdir(i)=k;
            xreg(i,:)=re;  i=i+1;
        elseif length(re)==3
            re = unique(re);
            re = sort(re);  xdir(i)=k;
            xreg(i,1)=re(1);  xreg(i,2)=re(2); i=i+1;
        end
        
    end % end if
end
xdir=xdir((xdir>0));
if isempty(xdir)
    xreg = 0;
else
    xreg = xreg(1:length(xdir),:);
end
return;
