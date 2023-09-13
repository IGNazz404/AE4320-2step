clc;
close all;
clear all;

syms x y z u v w phi theta psi wx wy wz

xgps = x;
ygps = y;
zgps = z;
ugps = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psi) - (v*cos(phi)-w*sin(phi))*sin(psi) + wx;
vgps = (u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psi) + (v*cos(phi)-w*sin(phi))*cos(psi) + wy;
wgps = -u*sin(theta) + (v*sin(phi)+w*cos(phi))*cos(theta) + wz;
phigps = phi;
thetagps = theta;
psigps = psi;
V = sqrt(u^2+v^2+w^2);
alfa = atan(w/u);
beta = atan(v/sqrt(u^2+w^2));

zobs = [xgps; ygps; zgps; ugps; vgps; wgps; phigps; thetagps; psigps; V; alfa; beta];

save dumpfile zobs
clear
load dumpfile
