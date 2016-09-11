% Purpose : Setup script, building operators, grid, metric, and connectivity tables.
% Definition of constants
Nfp = N+1; Np = (N+1)*(N+2)/2; Nfaces=3; 

% Compute nodal set
[x,y] = Nodes2D(N); % 等边三角形上的非等距插值点
[r,s] = xytors(x,y);% 标准三角形上的插值点

% Build reference element matrices 在标准单元上

V = Vandermonde2D(N,r,s); invV = inv(V);
MassMatrix = invV'*invV; % 局部质量矩阵
[Dr,Ds] = Dmatrices2D(N, r, s, V);
invMassMatrix = V*V'; % inverse of local matrix 

% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';% EToV:每个单元的三个节点的全局编号
x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));% 每一行指对应的每个影射单元上所有节点的横坐标
y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));% 每一行指对应的每个影射单元上所有节点的纵坐标

% find all the nodes that lie on each edge % 局部编号
fmask1   = find( abs(s+1) < NODETOL)'; % 第一条边上的节点
fmask2   = find( abs(r+s) < NODETOL)'; % 第二条边上的节点
fmask3   = find( abs(r+1) < NODETOL)'; % 第三条边上的节点
Fmask  = [fmask1;fmask2;fmask3]';
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);%所有边界点的横坐标
% Create surface integral terms
LIFT = Lift2D();

% calculate geometric factors
[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);

% calculate geometric factors
[nx, ny, sJ] = Normals2D();
Fscale = sJ./(J(Fmask,:));

gloinvmat = zeros(Np*K,Np*K); % inverse of global mass mat
for k=1:K
    gloinvmat((k-1)*Np+1:k*Np,(k-1)*Np+1:k*Np) = gloinvmat((k-1)*Np+1:k*Np,(k-1)*Np+1:k*Np) + (V*V'); 
end
% Build connectivity matrix
% EToE: EToE(k,i)为与单元k面i相联的单元整体编号
% EToF: EToF(k,i)为与单元k面i相联的局部面编号

[EToE, EToF] = tiConnect2D(EToV);
BCType = int8(not(EToE'-ones(Nfaces,1)*(1:K)))';
BCType=6*double(BCType);
% Build connectivity maps
BuildMaps2D;

% Compute weak operators (could be done in preprocessing to save time)
[Vr, Vs] = GradVandermonde2D(N, r, s);
Drw = (V*Vr')/(V*V'); Dsw = (V*Vs')/(V*V');
