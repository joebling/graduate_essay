function [VX, VY, K, EToV] = MeshGenDistMesh2D()
% function [VX, VY, K, EToV] = MeshGenDistMesh2D()
% Purpose : Generate 2D square mesh using DistMesh;
% By Allan P. Engsig-Karup
% Parameters to set/define
% fd Distance function for mesh boundary
% fh Weighting function for distributing elements
% h0 Characteristic length of elements
% Bbox Bounding box for mesh
% param Parameters to be used in function call with DistMesh
%%
% fd = inline('drectangle(p,-1,1,-1,1)','p');
% fh = @huniform;
% h0 = 0.5;
% Bbox = [-1 -1; 1 1];
% param = [];
% % Call distmesh
% [Vert,EToV]=distmeshnd(fd,fh,h0,Bbox,param);
% VX = Vert(:,1)'; VY = Vert(:,2)';
% Nv = length(VX); K = size(EToV,1);
% % Reorder elements to ensure counter clockwise orientation
% ax = VX(EToV(:,1)); ay = VY(EToV(:,1));
% bx = VX(EToV(:,2)); by = VY(EToV(:,2));
% cx = VX(EToV(:,3)); cy = VY(EToV(:,3));
% D = (ax-cx).*(by-cy)-(bx-cx).*(ay-cy);
% i = find(D<0);
% EToV(i,:) = EToV(i,[1 3 2]);
%%

str  = 'Rectangle';

% Characteristic length & Bounding box
h0   = 0.13;

% Domain parameters
xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;
Bbox = [xmin,ymin;xmax,ymax];
pfix = [xmin,ymin; xmin, ymax; xmax, ymin; xmax, ymax];

% Generating distance function
fd   = inline(['drectangle(p,' num2str(xmin) ',' num2str(xmax) ',' num2str(ymin) ',' num2str(ymax) ')'],'p');

% Element size function
choice = 4;
switch choice
    case 1 % south west corner refined
        fh   = inline(['min(0.5+sqrt((p(:,1)-' num2str(xmin) ').^2 + (p(:,2)-' num2str(ymin) ').^2 ), 1.4 )'],'p');
    case 2 % west side refined
        fh   = inline('4+2*p(:,1)','p');
    case 3 % middle refined
        fh   = inline(['min(0.5+sqrt((p(:,1)-' num2str((xmin+xmax)/2) ').^2 + (p(:,2)-' num2str((ymin+ymax)/2) ').^2 ), 1.4 )'],'p');
    case 4 % all corners refined
        fh   = inline(['min(2-1.1*sqrt((p(:,1)-' num2str((xmin+xmax)/2) ').^2 + (p(:,2)-' num2str((ymin+ymax)/2) ').^2 ), 1 )'],'p');
    case 5 % specific area refined
        fh   = inline('max(min(0.5*sqrt((p(:,1)-0.1).^2 + (p(:,2)-0.1).^2 ), 0.3 ),0.1)','p');
    case 6 % uniform mesh
        fh   = @huniform;
end

% Generate mesh:
disp('Generating mesh...')
[Vert,EToV]=distmesh2d(fd,fh,h0,Bbox,pfix);
disp('DONE!')

% Determine mesh parameters
VX   = Vert(:,1)';
VY   = Vert(:,2)';

% cd(CURRENTDIR)

EToV = Reorder(EToV,VX,VY);
Nv   = length(VX);
K    = size(EToV,1);

% Output mesh information:
disp('======================================================')
disp(['Mesh ....................... : ' str])
disp(['Number of elements ......... : ' num2str(K)])
disp(['Number of unique vertices .. : ' num2str(Nv)])



return