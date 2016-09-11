function rhsu = fdeRHS2D2(u,time)
Globals2D;

%% Define field differences at faces and impose Dirichlet BCs
% du = zeros(Nfp*Nfaces,K); du(:) = u(vmapM)-u(vmapP); %jump 
% du(mapD) = 2*(u(vmapD));
% % du(mapD) = 2*(u(vmapD) - exp(-time)*sin(pi*Fx(mapD)).*sin(pi*Fy(mapD))); 
% % du(mapD) = 2*(u(vmapD) - exp(-time)*sin(pi*(Fx(mapD)+1).^2).*sin(pi*(Fy(mapD)+1).^2)); 
% % Compute qx and qy, define differences and impose Neumann BC's
% [dudx,dudy] = Grad2D(u);
% 
% % Compute DG gradient with central fluxes
% fluxxu = nx.*du/2; qx = dudx - LIFT*(Fscale.*fluxxu); 
% fluxyu = ny.*du/2; qy = dudy - LIFT*(Fscale.*fluxyu);
% 
% % Compute minimum height of elements either side of each edge
% hmin = min(2*J(vmapP)./sJ(mapP), 2*J(vmapM)./sJ(mapM));
% tau = reshape(Np./hmin, Nfp*Nfaces, K); 
% 
% % Evaluate jumps in components of q at element interfaces
% dqx=zeros(Nfp*Nfaces,K); dqx(:)=(qx(vmapM)-qx(vmapP));  dqx(mapN) = 2*qx(vmapN);
% dqy=zeros(Nfp*Nfaces,K); dqy(:)=(qy(vmapM)-qy(vmapP));  dqy(mapN) = 2*qy(vmapN);
% 
% [dqxdx,dqxdy] = Grad2D(qx);dqxdx = dqxdx - LIFT*(Fscale.*(nx.*dqx)/2);
% [dqydx,dqydy] = Grad2D(qy);dqydy = dqydy - LIFT*(Fscale.*(ny.*dqy)/2);
% 
% % % Evaluate flux function
% % fluxq = (nx.*dqx + ny.*dqy)/2;
% % % fluxq = (nx.*dqx + ny.*dqy + tau.*du)/2;
% % 
% % % Compute right hand side
% % % divq = Div2D(qx,qy);
% % divq = divq - LIFT*(Fscale.*fluxq);
% % rhsu = reshape(globalFDX*dqxdx(:), Np, K) +  reshape(globalFDY*dqydy(:), Np, K)  + exp(-time)*f;
% % rhsu = reshape (gloinvmat*globalFIX*dqxdx(:),Np,K) + reshape(gloinvmat*globalFIY*dqydy(:), Np, K)  + exp(-time)*f;
% qqx = zeros(size(dqxdx));
% for k=1:length(xflag)
%     row = xflag(k,1);
%     col = xflag(k,2);
%     qqx(:,row) = qqx(:,row) + globalFIX(:,(k-1)*Np+1:k*Np)*dqxdx(:,col);
% end
% qqy = zeros(size(dqydy));
% for k=1:length(yflag)
%     row = yflag(k,1);
%     col = yflag(k,2);
%     qqy(:,row) = qqy(:,row) + globalFIY(:,(k-1)*Np+1:k*Np)*dqydy(:,col);
% end
% rhsu = invMassMatrix*qqx + invMassMatrix*qqy + exp(-time)*f;
% rhsu = reshape(globalFDX*(divq(:)), Np, K) + reshape(globalFDY*divq(:), Np, K) + exp(-time)*f;
%%
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapM)-u(vmapP); %jump 
% du(mapD) = 2*(u(vmapD));
du(mapD) = (u(vmapD) - exp(-time)*sin(pi*Fx(mapD)).*sin(pi*Fy(mapD))); 
% du(mapD) = 2*(u(vmapD) - exp(-time)*sin(pi*(Fx(mapD)+1).^2).*sin(pi*(Fy(mapD)+1).^2)); 
% Compute qx and qy, define differences and impose Neumann BC's
[dudx,dudy] = Grad2D(u);

% Compute DG gradient with central fluxes
fluxxu = nx.*du/2; qx = dudx - LIFT*(Fscale.*fluxxu); 
fluxyu = ny.*du/2; qy = dudy - LIFT*(Fscale.*fluxyu);

% Compute minimum height of elements either side of each edge
hmin = min(2*J(vmapP)./sJ(mapP), 2*J(vmapM)./sJ(mapM));
tau = reshape(Np./hmin, Nfp*Nfaces, K); 

qqx = zeros(size(qx));
for k=1:length(xflag)
    row = xflag(k,1);
    col = xflag(k,2);
    qqx(:,row) = qqx(:,row) + globalFIX(:,(k-1)*Np+1:k*Np)*qx(:,col);
end
qqy = zeros(size(qy));
for k=1:length(yflag)
    row = yflag(k,1);
    col = yflag(k,2);
    qqy(:,row) = qqy(:,row) + globalFIY(:,(k-1)*Np+1:k*Np)*qy(:,col);
end
qqx = invMassMatrix*qqx; qqy = invMassMatrix*qqy;

dqqx=zeros(Nfp*Nfaces,K); dqqx(:)=(qqx(vmapM)-qqx(vmapP));  dqqx(mapN) = 2*qqx(vmapN);
dqqy=zeros(Nfp*Nfaces,K); dqqy(:)=(qqy(vmapM)-qqy(vmapP));  dqqy(mapN) = 2*qqy(vmapN);

divq = Div2D(qqx,qqy);
fluxq = (nx.*dqqx + ny.*dqqy)/2;
rhsu = divq  - LIFT*(Fscale.*fluxq) + exp(-time)*f;

return;