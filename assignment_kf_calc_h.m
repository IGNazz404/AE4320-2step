 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   zpred = kf_calcHx(x) Calculates the output dynamics equation h(x,u,t) 
%
%   Author: Y.J.E. Prencipe
%   Student Number: 4777158
%   Course: AE4320 System Identification of Aerospace Vehicles
%   Place: Delft University of Technology, 2023
%   Email: y.j.e.prencipe@student.tudelft.nl
%   Version: 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zpred = assignment_kf_calc_h(x)
       
    % Output Equations (from assignmen eq. (5) and eq. (6))
    zpred = zeros(12, 1);

    % x_GPS; y_GPS; z_GPS
    zpred(1) = x(1); 
    zpred(2) = x(2); 
    zpred(3) = x(3);
    % u_GPS
    zpred(4) = (x(4)*cos(x(8)) + (x(5)*sin(x(7))+x(6)*cos(x(7)))*sin(x(8)))*cos(x(9)) - (x(5)*cos(x(7)) - x(6)*sin(x(7)))*sin(x(9)) + x(10);
    % v_GPS
    zpred(5) = (x(4)*cos(x(8)) + (x(5)*sin(x(7))+x(6)*cos(x(7)))*sin(x(8)))*sin(x(9)) + (x(5)*cos(x(7)) - x(6)*sin(x(7)))*cos(x(9)) + x(11);
    % w_GPS
    zpred(6) = -x(4)*sin(x(8)) + (x(5)*sin(x(7))+x(6)*cos(x(7)))*cos(x(8)) + x(12);
    % phi_GPS; theta_GPS; psi_GPS
    zpred(7) = x(7); 
    zpred(8) = x(8); 
    zpred(9) = x(9);
    % V; alfa; beta
    zpred(10) = sqrt(x(4)^2+x(5)^2+x(6)^2); 
    zpred(11) = atan(x(6)/x(4)); 
    zpred(12) = atan(x(5)/sqrt(x(4)^2+x(6)^2));

end
