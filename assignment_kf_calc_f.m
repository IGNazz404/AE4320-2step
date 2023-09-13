%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   xdot = kf_calcFx(x) Calculates the system dynamics equation f(x,u,t) 
%
%   Author: Y.J.E. Prencipe
%   Student Number: 4777158
%   Course: AE4320 System Identification of Aerospace Vehicles
%   Place: Delft University of Technology, 2023
%   Email: y.j.e.prencipe@student.tudelft.nl
%   Version: 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = assignment_kf_calc_f(t, x, u)

    xdot = zeros(18, 1);
    g = 9.7897; %gravitywgs84(7500, 52+23/60);
    % g = 9.7837; % proper value, not used for reporting, found using the
    % state prediction error --> aircraft flies at 9500m

    % % system dynamics go here!
    xdot(1) = (x(4)*cos(x(8)) + (x(5)*sin(x(7)) + x(6)*cos(x(7))) * sin(x(8)))*cos(x(9)) - (x(5)*cos(x(7)) - x(6)*sin(x(7)))*sin(x(9)) + x(10);
    xdot(2) = (x(4)*cos(x(8)) + (x(5)*sin(x(7)) + x(6)*cos(x(7))) * sin(x(8)))*sin(x(9)) + (x(5)*cos(x(7)) - x(6)*sin(x(7)))*cos(x(9)) + x(11);
    xdot(3) = -x(4)*sin(x(8)) + (x(5)*sin(x(7)) + x(6)*cos(x(7)))*cos(x(8)) + x(12);
    xdot(4) = u(1)-x(13) - g*sin(x(8)) + (u(6)-x(18))*x(5) - (u(5)-x(17))*x(6);
    xdot(5) = u(2)-x(14) + g*cos(x(8))*sin(x(7)) + (u(4)-x(16))*x(6) - (u(6)-x(18))*x(4);
    xdot(6) = u(3)-x(15) + g*cos(x(8))*cos(x(7)) + (u(5)-x(17))*x(4) - (u(4)-x(16))*x(5);
    xdot(7) = u(4)-x(16) + (u(5)-x(17))*sin(x(7))*tan(x(8)) + (u(6)-x(18))*cos(7)*tan(x(8));
    xdot(8) = (u(5)-x(17))*cos(x(7)) - (u(6)-x(18))*sin(x(7));
    xdot(9) = (u(5)-x(17))*sin(x(7))/cos(x(8)) + (u(6)-x(18))*cos(x(7))/cos(x(8));
    xdot(10) = 0;
    xdot(11) = 0;
    xdot(12) = 0;   % Derivatives of the wind speed are 0 since the wind speed components are assumed to be constants
    xdot(13) = 0;
    xdot(14) = 0;
    xdot(15) = 0;
    xdot(16) = 0;
    xdot(17) = 0;
    xdot(18) = 0;   % Derivatives of input biases is 0 since theyre constants

end
