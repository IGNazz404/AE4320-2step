%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Structure of the Alternative Model for Least Squares Estimation
%
%   Author: Y.J.E. Prencipe, based on C.C. de Visser
%   Student Number: 4777158
%   Course: AE4320 System Identification of Aerospace Vehicles
%   Place: Delft University of Technology, 2023
%   Email: y.j.e.prencipe@student.tudelft.nl
%   Version: 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pdotnorm = (pdot./(2*Vm))'.*b; pdotnorm = (pdotnorm./(2*Vm')).*b;
qdotnorm = (qdot./(2*Vm))'.*c; qdotnorm = (qdotnorm./(2*Vm')).*c;
rdotnorm = (rdot./(2*Vm))'.*b; rdotnorm = (rdotnorm./(2*Vm')).*b;
Axc = (Axc./(2*Vm))'.*b;
Ayc = (Ayc./(2*Vm))'.*c;
Azc = (Azc./(2*Vm))'.*b;

M2.Cx = [ones(N,1), alpham, alpham.^2, qnorm, de, Axc, u, qdotnorm];
M2.Cy = [ones(N,1), betam, pnorm, rnorm, da, dr, Ayc, v, pdotnorm, rdotnorm];
M2.Cz = [ones(N,1), alpham, qnorm, de, Azc, u, w, qdotnorm];
M2.Cl = [ones(N,1), betam, pnorm, rnorm, da, dr, pdotnorm, rdotnorm];
M2.Cm = [ones(N,1), alpham, qnorm, de, qdotnorm];
M2.Cn = [ones(N,1), betam, pnorm, rnorm, da, dr, pdotnorm, rdotnorm];
