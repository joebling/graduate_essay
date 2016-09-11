function rhsu = fdeLDGHS2D(u,time)
Globals2D;

%% DG for 2D diffusion problem
% Define field differences at faces and impose Dirichlet BCs
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapM)-u(vmapP); %jump 
% du(mapD) = 2*(u(vmapD));

% impose boundary condition -- Dirichlet BC¡¯s
du(mapD) = (u(vmapD) - exp(-time)*((Fx(mapD).^2 - 1).^3).*((Fy(mapD).^2 - 1).^3));

% Compute qx and qy, define differences and impose Neumann BC's
[dudx,dudy] = Grad2D(u);

% Compute DG gradient with central fluxes
fluxxu = nx.*(du + nx.*du)/2; 
fluxyu = ny.*(du + ny.*du)/2;
% fluxxu = (1+(abs(nx)>1e-14)).*du/ 2; qx = dudx - LIFT*(Fscale.*(nx.*fluxxu)); 
% fluxyu = (1+(abs(ny)>1e-14)).*du/ 2; qy = dudy - LIFT*(Fscale.*(ny.*fluxyu));


% % compute qx and qy
qx = dudx - LIFT*(Fscale.*fluxxu); 
qy = dudy - LIFT*(Fscale.*fluxyu);

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
% qqx = invMassMatrix*qqx; 
% qqy = invMassMatrix*qqy;
% tau = 4*N*(N+1);
% qqx = invMassMatrix*qqx - LIFT*(Fscale.*(tau.*nx.*du)); 
% qqy = invMassMatrix*qqy - LIFT*(Fscale.*(tau.*ny.*du));

dqqx=zeros(Nfp*Nfaces,K); dqqx(:)=(qqx(vmapM)-qqx(vmapP));  dqqx(mapD) = 0;
dqqy=zeros(Nfp*Nfaces,K); dqqy(:)=(qqy(vmapM)-qqy(vmapP));  dqqy(mapD) = 0;

divq = Div2D(qqx,qqy);
fluxq = (nx.*(1 - nx).*dqqx + ny.*(1 - ny).*dqqy)/2;
% fluxq = (nx.*(1 - (abs(nx)>1e-14)).*dqqx + ny.*(1 - (abs(ny)>1e-14)).*dqqy)/2;
rhsu = divq  - LIFT*(Fscale.*fluxq) + exp(-time)*f;

%% DG for 2D diffusion problem
% Define field differences at faces and impose Dirichlet BCs
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
% - LIFT*(Fscale.*(nx.*(ddqxdx + nx.*ddqxdx/2)/2 )) - LIFT*(Fscale.*(ny.*(ddqydy + ny.*ddqydy/2)/2 ))  
return;
