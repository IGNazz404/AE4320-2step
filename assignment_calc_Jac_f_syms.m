%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Jacobian of function f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

syms x y z u v w phi theta psi wx wy wz lx ly lz lp lq lr Ax Ay Az p q r g


dx =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi)-(v*cos(phi)-w*sin(phi))*sin(psi)+wx;
dy =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi)+(v*cos(phi)-w*sin(phi))*cos(psi)+wy;
dz = -u*sin(theta) +(v*sin(phi)+w*cos(phi))*cos(theta)+wz;
du = Ax-lx - g*sin(theta)+(r-lr)*v - (q-lq)*w;
dv = Ay-ly + g*cos(theta)*sin(phi)+(p-lp)*w-(r-lr)*u;
dw = Az-lz + g*cos(theta)*cos(phi)+(q-lq)*u-(p-lp)*v;
dphi = (p-lp) + (q-lq)*sin(phi)*tan(theta)+(r-lr)*cos(phi)*tan(theta);
dtheta = (q-lq)*cos(phi)-(r-lr)*sin(phi);
dpsi = (q-lq)*sin(phi)/cos(theta)+(r-lr)*cos(phi)/cos(theta);
dW = 0;
dBias = 0;

dx = [dx;dy;dz;du;dv;dw;dphi;dtheta;dpsi;dW;dW;dW;dBias;dBias;dBias;dBias;dBias;dBias];

DFx = jacobian(dx, [x y z u v w phi theta psi wx wy wz lx ly lz lp lq lr]);
