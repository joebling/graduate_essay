% discontinous Galerkin method for solvng FDE on domain with non-uniform
% meshes, the exact solution is 
% u(x,y,t) = exp(-t)*(x^2 - 1)^3*(y^2-1)^3
clear;

tic;
Globals2D;

xorder = 1.1;
yorder = 1.1;


% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D('Maxwell025.neu');
[Nv, VX, VY, K, EToV] = equaltriangle(3);
% load mesh97;
display(K);

% for N=1:3
N = 1;
finalresult = 0;
% Polynomial order used for approximation 
% N = 2;

% initial point
Ia = -1; Ic = -1;
 

% Initialize solvr and construct grid and metric
StartUp2D;

% set up boundary conditions
BuildBCMaps2D;

% initial condition
u = ((x.^2-1).^3).*((y.^2-1).^3);

Finaltime = 0.01;
u = fde2D(u, Finaltime);

u = reshape(u, Np, K); % numerical solution

eu = exp(-Finaltime)*((x.^2-1).^3).*((y.^2-1).^3);
maxerror = max(max(abs(u-eu)));
display(maxerror);

L2error = 0;
for k=1:K
    B = [VX(EToV(k,2)) - VX(EToV(k,1)),VX(EToV(k,3)) - VX(EToV(k,1));
         VY(EToV(k,2)) - VY(EToV(k,1)),VY(EToV(k,3)) - VY(EToV(k,1))];
     detB = det(B);
     for i =1:length(trigauss)
         % trigauss: Guass integral points and weights on standard triangle
         % saved in Globals2D;
         eta = B*(trigauss(i,1:2))' + [VX(EToV(k,1));VY(EToV(k,1))];
         L2error = L2error + detB*( u(:,k)'*LagrangeInt2DTRI(eta(1),eta(2),k) - exp(-Finaltime)*(((eta(1)).^2 - 1).^3).*(((eta(2)).^2 - 1).^3) )^2*trigauss(i,3);
     end
end
L2error = sqrt(L2error);
display(L2error);
% end

% figure; PlotField2D(round((N-1)*10), x, y, u); drawnow; 
toc;
