function re=LagrangeInt2DTRI(x1,y1,m)
% Lagrange interpolation polynomial on standard triangle region, (-1,-1),(1,-1),(-1,1).
% In order to use this function, we must call for Vandermonde2D matrix.
% point(x1,y1) belongs to original triangle in physical space.
% m is the index of element affected by collocation point.
Globals2D;
re=zeros(Np,1);
sk=1;
vx = VX(EToV(m,:));  vy = VY(EToV(m,:));
[r1,s1] = xy2rs(x1,y1,vx,vy); %transform (x1,y1) to standard triangle.
[a, b] = rstoab(r1, s1); %transform (r1,s1) to unit sqaure .
for i=0:N
  for j=0:N - i
    re(sk) = Simplex2DP(a,b,i,j);
    sk = sk+1;
  end
end
re=inv(V')*re;
return;
