function [u] = fde2D(u, FinalTime)

Globals2D;
time = 0;

rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = (min(dtscale)*rmin)^2*0.5;
dt = 1.0e-5;
f = sourcefun();
% [globalFIX,globalFIY] = globalFIM(); % form the fractional stiffness matrix 
[globalFIX,xflag,globalFIY,yflag] = globalFIM2(); % low-store way
% load N2K18
% load N4M05;
% load data2;
% load mat1;
resu = zeros(Np,K);
while (time<FinalTime)
  if(time+dt>FinalTime), dt = FinalTime-time; end
%   u = u + dt*fdeRHS2D(u,time);
   for intrk = 1:5
      timelocal = time +rk4c(intrk)*dt;
       % initiate and increment Runge-Kutta residuals
      rhsu = fdeLDGHS2D(u,timelocal);
%       rhsu = fdeLDGHSPenalty2D(u,timelocal);

%       rhsu = fdeRHS2D(u,timelocal);
      resu = rk4a(intrk)*resu + dt*rhsu;
      u = u + rk4b(intrk)*resu;
   end
   % Increment time
   time = time+dt;
end
% u(vmapD) = exp(-FinalTime)*((Fx(mapD).^2 - 1).^3).*((Fy(mapD).^2 - 1).^3);
return;
