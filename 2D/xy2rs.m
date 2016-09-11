function [r1,s1] = xy2rs(x1,y1,vx,vy)
% Globals2D;
% Purpose : From (x1,y1) in physical triangle to (r,s) coordinates in standard triangle

lemda1 = ((vx(1)-vx(2))*(y1-vy(2)) - (vy(1)-vy(2))*(x1-vx(2)))/((vx(1)-vx(2))*(vy(3)-vy(2)) - (vy(1)-vy(2))*(vx(3)-vx(2)));
lemda2 = ((vx(3)-vx(2))*(y1-vy(2)) - (vy(3)-vy(2))*(x1-vx(2)))/((vx(3)-vx(2))*(vy(1)-vy(2)) - (vy(3)-vy(2))*(vx(1)-vx(2)));
lemda3 = 1 - lemda1 - lemda2;

r1 = -lemda2 + lemda3 - lemda1;
s1 = -lemda2 - lemda3 + lemda1;
return;
