%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Find Jacobian of measurement equation
%
%   Author: Y.J.E. Prencipe
%   Student Number: 4777158
%   Course: AE4320 System Identification of Aerospace Vehicles
%   Place: Delft University of Technology, 2023
%   Email: y.j.e.prencipe@student.tudelft.nl
%   Version: 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
clear all;

observation_eqs;
syms x y z u v w phi theta psi wx wy wz

DHxgps = [diff(zobs(1),x), diff(zobs(1),y), diff(zobs(1),z), diff(zobs(1),u), diff(zobs(1),v), diff(zobs(1),w), diff(zobs(1),phi), diff(zobs(1),theta), diff(zobs(1),psi), ...
    diff(zobs(1),wx), diff(zobs(1),wy), diff(zobs(1),wz)];

DHygps = [diff(zobs(2),x), diff(zobs(2),y), diff(zobs(2),z), diff(zobs(2),u), diff(zobs(2),v), diff(zobs(2),w), diff(zobs(2),phi), diff(zobs(2),theta), diff(zobs(2),psi), ...
    diff(zobs(2),wx), diff(zobs(2),wy), diff(zobs(2),wz)];

DHzgps = [diff(zobs(3),x), diff(zobs(3),y), diff(zobs(3),z), diff(zobs(3),u), diff(zobs(3),v), diff(zobs(3),w), diff(zobs(3),phi), diff(zobs(3),theta), diff(zobs(3),psi), ...
    diff(zobs(3),wx), diff(zobs(3),wy), diff(zobs(3),wz)];

DHugps = [diff(zobs(4),x), diff(zobs(4),y), diff(zobs(4),z), diff(zobs(4),u), diff(zobs(4),v), diff(zobs(4),w), diff(zobs(4),phi), diff(zobs(4),theta), diff(zobs(4),psi), ...
    diff(zobs(4),wx), diff(zobs(4),wy), diff(zobs(4),wz)];

DHvgps = [diff(zobs(5),x), diff(zobs(5),y), diff(zobs(5),z), diff(zobs(5),u), diff(zobs(5),v), diff(zobs(5),w), diff(zobs(5),phi), diff(zobs(5),theta), diff(zobs(5),psi), ...
    diff(zobs(5),wx), diff(zobs(5),wy), diff(zobs(5),wz)];

DHwgps = [diff(zobs(6),x), diff(zobs(6),y), diff(zobs(6),z), diff(zobs(6),u), diff(zobs(6),v), diff(zobs(6),w), diff(zobs(6),phi), diff(zobs(6),theta), diff(zobs(6),psi), ...
    diff(zobs(6),wx), diff(zobs(6),wy), diff(zobs(6),wz)];

DHphigps = [diff(zobs(7),x), diff(zobs(7),y), diff(zobs(7),z), diff(zobs(7),u), diff(zobs(7),v), diff(zobs(7),w), diff(zobs(7),phi), diff(zobs(7),theta), diff(zobs(7),psi), ...
    diff(zobs(7),wx), diff(zobs(7),wy), diff(zobs(7),wz)];

DHthetagps = [diff(zobs(8),x), diff(zobs(8),y), diff(zobs(8),z), diff(zobs(8),u), diff(zobs(8),v), diff(zobs(8),w), diff(zobs(8),phi), diff(zobs(8),theta), diff(zobs(8),psi), ...
    diff(zobs(8),wx), diff(zobs(8),wy), diff(zobs(8),wz)];

DHpsigps = [diff(zobs(9),x), diff(zobs(9),y), diff(zobs(9),z), diff(zobs(9),u), diff(zobs(9),v), diff(zobs(9),w), diff(zobs(9),phi), diff(zobs(9),theta), diff(zobs(9),psi), ...
    diff(zobs(9),wx), diff(zobs(9),wy), diff(zobs(9),wz)];

DHV = [diff(zobs(10),x), diff(zobs(10),y), diff(zobs(10),z), diff(zobs(10),u), diff(zobs(10),v), diff(zobs(10),w), diff(zobs(10),phi), diff(zobs(10),theta), diff(zobs(10),psi), ...
    diff(zobs(10),wx), diff(zobs(10),wy), diff(zobs(10),wz)];

DHalfa = [diff(zobs(11),x), diff(zobs(11),y), diff(zobs(11),z), diff(zobs(11),u), diff(zobs(11),v), diff(zobs(11),w), diff(zobs(11),phi), diff(zobs(11),theta), diff(zobs(11),psi), ...
    diff(zobs(11),wx), diff(zobs(11),wy), diff(zobs(11),wz)];

DHbeta = [diff(zobs(12),x), diff(zobs(12),y), diff(zobs(12),z), diff(zobs(12),u), diff(zobs(12),v), diff(zobs(12),w), diff(zobs(12),phi), diff(zobs(12),theta), diff(zobs(12),psi), ...
    diff(zobs(12),wx), diff(zobs(12),wy), diff(zobs(12),wz)];


Dhx = [DHxgps; DHygps; DHzgps; DHugps; DHvgps; DHwgps; DHphigps; DHthetagps; DHpsigps; DHV; DHalfa; DHbeta];

save dumpfile Dhx
clear
load dumpfile
