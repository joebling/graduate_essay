 function rhsu = fdeLDGHSPenalty2D(u,time)
Globals2D;

%% DG for 2D diffusion problem
% Define field differences at faces and impose Dirichlet BCs
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapM)-u(vmapP); %jump 
% du(mapD) = 2*(u(vmapD));
du(mapD) = (u(vmapD) - exp(-time)*((Fx(mapD).^2 - 1).^3).*((Fy(mapD).^2 - 1).^3)); 


% Compute qx and qy, define differences and impose Neumann BC's
[dudx,dudy] = Grad2D(u);

% Compute DG gradient with central fluxes
fluxxu = nx.*(du + nx.*du)/2; qx = dudx - LIFT*(Fscale.*fluxxu); 
fluxyu = ny.*(du + ny.*du)/2; qy = dudy - LIFT*(Fscale.*fluxyu);
% fluxxu = (1+nx).*du/ 2; qx = dudx - LIFT*(Fscale.*(nx.*fluxxu)); 
% fluxyu = (1+ny).*du/ 2; qy = dudy - LIFT*(Fscale.*(ny.*fluxyu));

dqx=zeros(Nfp*Nfaces,K); dqx(:)=(qx(vmapM)-qx(vmapP));  dqx(mapN) = 2*qx(vmapN);
dqy=zeros(Nfp*Nfaces,K); dqy(:)=(qy(vmapM)-qy(vmapP));  dqy(mapN) = 2*qy(vmapN);

hmin = min(2*J(vmapP)./sJ(mapP), 2*J(vmapM)./sJ(mapM));
taux = 0.1*reshape((N*(N+1))^2*xorder./hmin, Nfp*Nfaces, K); 
tauy = 0.1*reshape((N*(N+1))^2*yorder./hmin, Nfp*Nfaces, K); 
% tau  = reshape(Np./hmin, Nfp*Nfaces, K);  
tau = 1;
% tau = (N*(N+1))^2/(min(hmin)/2);
% tau = 0;

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

qqx = invMassMatrix*qqx - LIFT*(Fscale.*(tau.*nx.*dqx)); 
qqy = invMassMatrix*qqy - LIFT*(Fscale.*(tau.*ny.*dqy));
% qqx(vmapD) = qqx(vmapD) - tau.*du(mapD);
% qqy(vmapD) = qqy(vmapD) - tau.*du(mapD);

% qqx = invMassMatrix*qqx;
% qqy = invMassMatrix*qqy;


dqqx=zeros(Nfp*Nfaces,K); dqqx(:)=(qqx(vmapM)-qqx(vmapP));  dqqx(mapN) = 2*qqx(vmapN);
dqqy=zeros(Nfp*Nfaces,K); dqqy(:)=(qqy(vmapM)-qqy(vmapP));  dqqy(mapN) = 2*qqy(vmapN);

divq = Div2D(qqx,qqy);
fluxq = (nx.*(1 - nx).*dqqx + ny.*(1 - ny).*dqqy)/2;
% rhsu = divq  - LIFT*(Fscale.*(fluxq+tau.*du)) + exp(-time)*f;
rhsu = divq  - LIFT*(Fscale.*(fluxq)) + exp(-time)*f;

%% DG for 2D diffusion problem
% % Define field differences at faces and impose Dirichlet BCs
% du = zeros(Nfp*Nfaces,K); du(:) = u(vmapM)-u(vmapP); %jump 
% du(mapD) = 2*(u(vmapD));
% % du(mapD) = 2*(u(vmapD) - exp(-time)*sin(pi*Fx(mapD)).*sin(pi*Fy(mapD))); 
% 
% % Compute qx and qy, define differences and impose Neumann BC's
% [dudx,dudy] = Grad2D(u);
% 
% % Compute DG gradient with central fluxes
% fluxxu = nx.*(du + nx.*du)/2; qx = dudx - LIFT*(Fscale.*fluxxu); 
% fluxyu = ny.*(du + ny.*du)/2; qy = dudy - LIFT*(Fscale.*fluxyu);
% % fluxxu = (1+nx).*du/ 2; qx = dudx - LIFT*(Fscale.*(nx.*fluxxu)); 
% % fluxyu = (1+ny).*du/ 2; qy = dudy - LIFT*(Fscale.*(ny.*fluxyu));
% 
% % Evaluate jumps in components of q at element interfaces
% dqx=zeros(Nfp*Nfaces,K); dqx(:)=(qx(vmapM)-qx(vmapP));  dqx(mapN) = 2*qx(vmapN);
% dqy=zeros(Nfp*Nfaces,K); dqy(:)=(qy(vmapM)-qy(vmapP));  dqy(mapN) = 2*qy(vmapN);
% 
% % Evaluate flux function
% % fluxq = (nx.*(dqx - nx.*dqx/2) + ny.*(dqy - ny.*dqy/2))/2;
% % fluxqx = (1 - nx).*dqx/ 2;
% % fluxqy = (1 - ny).*dqy/ 2;
% % Compute right hand side
% [dqxdx,dqxdy] = Grad2D(qx); dqxdx = dqxdx - LIFT*(Fscale.*(nx.*(dqx - nx.*dqx))/2);
% [dqydx,dqydy] = Grad2D(qy); dqydy = dqydy - LIFT*(Fscale.*(ny.*(dqy - ny.*dqy))/2);
% 
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
% rhsu = reshape(globalFDX*divq(:), Np, K) +  reshape(globalFDY*divq(:), Np, K)  + exp(-time)*f;
% %- LIFT*(Fscale.*(nx.*(ddqxdx + nx.*ddqxdx/2)/2 )) - LIFT*(Fscale.*(ny.*(ddqydy + ny.*ddqydy/2)/2 ))  
return;
