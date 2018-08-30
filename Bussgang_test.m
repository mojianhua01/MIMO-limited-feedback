clc;
close all;
[x ,y] = meshgrid(-10:0.1:10, -10:0.1:10);
r = 0.9;
z = 1/2/pi/sqrt(1-r^2) * exp( - (x.^2 + y.^2 - 2 *r.*x.*y)./(2*(1-r^2)) );

figure,
surf(x, y, z)

figure, 
contourf(x,y,z)

clear;
syms r x y;
int(1/2/pi/sqrt(1-r^2) * exp( - (x^2 + y^2 - 2 *r*x*y)/(2*(1-r^2)) ), x)
