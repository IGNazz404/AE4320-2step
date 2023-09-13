%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Derivatives of EOM
%
%   Author: Y.J.E. Prencipe
%   Student Number: 4777158
%   Course: AE4320 System Identification of Aerospace Vehicles
%   Place: Delft University of Technology, 2023
%   Email: y.j.e.prencipe@student.tudelft.nl
%   Version: 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms x y z u v w phi theta psi wx wy wz lx ly lz lp lq lr ax ay az p q r g

xdot = (u*cos(theta) + (v*sin(phi) + w*cos(phi)))*cos(psi) - (v*cos(phi) - w*sin(phi))*sin(psi) + wx;
ydot = (u*cos(theta) + (v*sin(phi) + w*cos(phi)))*sin(psi) + (v*cos(phi) - w*sin(phi))*cos(psi) + wy;
zdot = -u*sin(theta) + (v*sin(phi) + w*cos(phi))*cos(psi) + wz;
udot = (ax + lx) - g*sin(theta) + (r+lr)*v - (q+lq)*w;
vdot = (ay + ly) + g*cos(theta)*sin(phi) + (p+lp)*w - (r+lr)*u;
wdot = (az + lz) + g*cos(theta)*cos(phi) + (q+lq)*u - (p+lp)*v;
phidot = (p+lp) + (q+lq)*sin(phi)*tan(theta) + (r+lr)*cos(phi)*tan(theta);
thetadot = (q+lq)*cos(phi) - (r+lr)*sin(phi);
psidot = (q+lq+wq)*sin(phi)/cos(theta) + (r+lr+wr)*cos(phi)/cos(theta);
wxdot = 0;
wydot = 0;
wzdot = 0;
lxdot = 0;
lydot = 0;
lzdot = 0;
lpdot = 0;
lqdot = 0;
lrdot = 0;

xdot = [xdot, ydot, zdot, udot, vdot, wdot, phidot, thetadot, psidot, wxdot, wydot, wzdot, lxdot, lydot, lzdot, lpdot, lqdot, lrdot];

save dumpfile xdot
clear
load dumpfile
