%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   R^2 Computation Function
%
%   Author: Y.J.E. Prencipe
%   Student Number: 4777158
%   Course: AE4320 System Identification of Aerospace Vehicles
%   Place: Delft University of Technology, 2023
%   Email: y.j.e.prencipe@student.tudelft.nl
%   Version: 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rsquared = my_Rsquared_coeff(y,pred)

    %R2 correlation coefficient computation
    Rsquared = 1 - sum((y - pred).^2)/sum((y - mean(y)).^2);

end