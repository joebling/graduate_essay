% Purpose : Setup script, building operators, grid, metric, and connectivity tables.
% Definition of constants
Nfp = N+1; Np = (N+1)*(N+2)/2; Nfaces=3; 

% Compute nodal set
[x,y] = Nodes2D(N); % �ȱ��������ϵķǵȾ��ֵ��
[r,s] = xytors(x,y);% ��׼�������ϵĲ�ֵ��

% Build reference element matrices �ڱ�׼��Ԫ��

V = Vandermonde2D(N,r,s); invV = inv(V);
MassMatrix = invV'*invV; % �ֲ���������
[Dr,Ds] = Dmatrices2D(N, r, s, V);
invMassMatrix = V*V'; % inverse of local matrix 

% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';% EToV:ÿ����Ԫ�������ڵ��ȫ�ֱ��
x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));% ÿһ��ָ��Ӧ��ÿ��Ӱ�䵥Ԫ�����нڵ�ĺ�����
y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));% ÿһ��ָ��Ӧ��ÿ��Ӱ�䵥Ԫ�����нڵ��������

% find all the nodes that lie on each edge % �ֲ����
fmask1   = find( abs(s+1) < NODETOL)'; % ��һ�����ϵĽڵ�
fmask2   = find( abs(r+s) < NODETOL)'; % �ڶ������ϵĽڵ�
fmask3   = find( abs(r+1) < NODETOL)'; % ���������ϵĽڵ�
Fmask  = [fmask1;fmask2;fmask3]';
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);%���б߽��ĺ�����
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
% EToE: EToE(k,i)Ϊ�뵥Ԫk��i�����ĵ�Ԫ������
% EToF: EToF(k,i)Ϊ�뵥Ԫk��i�����ľֲ�����

[EToE, EToF] = tiConnect2D(EToV);
BCType = int8(not(EToE'-ones(Nfaces,1)*(1:K)))';
BCType=6*double(BCType);
% Build connectivity maps
BuildMaps2D;

% Compute weak operators (could be done in preprocessing to save time)
[Vr, Vs] = GradVandermonde2D(N, r, s);
Drw = (V*Vr')/(V*V'); Dsw = (V*Vs')/(V*V');
