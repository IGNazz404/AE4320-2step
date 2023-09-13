%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Jacobian of function h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

syms x y z u v w phi theta psi wx wy wz lx ly lz lp lq lr


xm = x;
ym = y;
zm = z;
um = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi)-(v*cos(phi)-w*sin(phi))*sin(psi)+wx;
vm = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi)+(v*cos(phi)-w*sin(phi))*cos(psi)+wy;
wm = -u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+wz;
phim = phi;
thetam = theta;
psim = psi;
vtasm = sqrt(u^2+v^2+w^2);
alpham = atan(w/u);
betam = atan(v/sqrt(u^2+w^2));

zm = [xm;ym;zm;um;vm;wm;phim;thetam;psim;vtasm;alpham;betam];

Hx = jacobian(zm, [x y z u v w phi theta psi wx wy wz lx ly lz lp lq lr]);
