% Driver script for solving the 2D vacuum Maxwell's equations on TM form
Globals2D;

% Polynomial order used for approximation 
N = 3;

% Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Maxwell05.neu');

% Initialize solver and construct grid and metric
StartUp2D;

% Set initial conditions
mmode = 1; nmode = 1;
Ez = sin(mmode*pi*x).*sin(nmode*pi*y); Hx = zeros(Np, K); Hy = zeros(Np, K);

% Solve Problem
FinalTime = 1;
[Hx,Hy,Ez,time] = Maxwell2D(Hx,Hy,Ez,FinalTime);

omega = pi*sqrt(mmode^2+nmode^2);
% exact solution
hx = -pi*nmode/omega*sin(mmode*pi*x).*cos(nmode*pi*y)*sin(omega*FinalTime);
hy = pi*mmode/omega*cos(mmode*pi*x).*sin(nmode*pi*y)*sin(omega*FinalTime);
ez = sin(mmode*pi*x).*sin(nmode*pi*y)*cos(omega*FinalTime);

maxerror1 = max(max(abs(Hx-hx)));
maxerror2 = max(max(abs(Hy-hy)));
maxerror3 = max(max(abs(Ez-ez)));

display(maxerror1);
display(maxerror2);
display(maxerror3);
