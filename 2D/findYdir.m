function [ydir,yreg]=findYdir(px,py)
Globals2D;
ydir=zeros(K,1);
yreg=Inf(K,2);
i=1;
for k=1:K
    minx = min( VX(EToV(k,:))); maxx = max( VX(EToV(k,:)));
    if px>=minx && px<=maxx
        re= [];
        vx = VX(EToV(k,:));  vy = VY(EToV(k,:));
        if (vx(1) - px)*(vx(2) - px)<=0
%             if abs(vx(2)-vx(1))<=NODETOL
%                 re = [re,py];
%             else
                re = [re,(vy(2) - vy(1))/(vx(2)-vx(1))*(px-vx(1)) + vy(1)];
%             end
        end
        if (vx(3) - px)*(vx(2) - px)<=0
%             if abs(vx(2)-vx(3))<NODETOL
%                 re = [re,py];
%             else
                re = [re,(vy(2) - vy(3))/(vx(2)-vx(3))*(px-vx(3)) + vy(3)];
%             end
        end
        if (vx(1) - px)*(vx(3) - px)<=0
%             if abs(vx(3)-vx(1))<NODETOL
%                 re = [re,py];
%             else
                re = [re,(vy(3) - vy(1))/(vx(3)-vx(1))*(px-vx(1)) + vy(1)];
%             end
        end
        re = re(re<=py+NODETOL);
       %%
        if length(re)==1
            ydir(i)=k;
            yreg(i,1) = re; i=i+1; 
        elseif length(re)==2
            re = sort(re);  ydir(i)=k;
            yreg(i,1) = re(1); yreg(i,2) = re(2); i=i+1;
        elseif length(re)==3
            re = unique(re);re = sort(re);ydir(i)=k;
            yreg(i,1) = re(1); yreg(i,2) = re(2); i=i+1;
        end
    end
end
ydir=ydir((ydir>0));
if isempty(ydir)
    yreg = 0;
else
    yreg = yreg(1:length(ydir),:);
end
return;
