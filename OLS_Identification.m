%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OLS Identification Model
%
%   Author: Y.J.E. Prencipe, based on C.C. de Visser
%   Student Number: 4777158
%   Course: AE4320 System Identification of Aerospace Vehicles
%   Place: Delft University of Technology, 2023
%   Email: y.j.e.prencipe@student.tudelft.nl
%   Version: 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Concatenated X Matrices for Original Model Structure
OLS.ID.X.XCx_datasets = [];
OLS.ID.X.XCy_datasets = [];
OLS.ID.X.XCz_datasets = [];
OLS.ID.X.XCl_datasets = [];
OLS.ID.X.XCm_datasets = [];
OLS.ID.X.XCn_datasets = [];
% Define Concatenated X Matrices for Alternative Model Structure
OLS.ID.M2.X.XCx = [];
OLS.ID.M2.X.XCy = [];
OLS.ID.M2.X.XCz = [];
OLS.ID.M2.X.XCl = [];
OLS.ID.M2.X.XCm = [];
OLS.ID.M2.X.XCn = [];
% Define Concatenated Y Matrices
OLS.ID.Y.YCx_datasets = [];
OLS.ID.Y.YCy_datasets = [];
OLS.ID.Y.YCz_datasets = [];
OLS.ID.Y.YCl_datasets = [];
OLS.ID.Y.YCm_datasets = [];
OLS.ID.Y.YCn_datasets = [];

for i=1:4

    % Run IEKF for Selected File
    fileNumber = i;
    dataFile = ["da3211_2.mat", "dadoublet_1.mat", "de3211_1.mat", "dr3211_1.mat", "dedoublet_1.mat", "dr3211_2.mat"];
    dataSet  = dataFile(fileNumber);
    load(dataSet);
    assignmentMultiStateIteratedEKF;
    
    % Bias Estimates
    lambdax = mean(XX_k1_k1(13,:), "all"); 
    lambday = mean(XX_k1_k1(14,:), "all"); 
    lambdaz = mean(XX_k1_k1(15,:), "all"); 
    lambdap = mean(XX_k1_k1(16,:), "all"); 
    lambdaq = mean(XX_k1_k1(17,:), "all"); 
    lambdar = mean(XX_k1_k1(18,:), "all"); 
    
    % Measured Inputs U_k Corrected for Bias (hence they still include Noise)
    Axc = U_k(1,:) - lambdax;   % c denotes 'corrected'
    Ayc = U_k(2,:) - lambday;
    Azc = U_k(3,:) - lambdaz;
    pc  = U_k(4,:) - lambdap;
    qc  = U_k(5,:) - lambdaq;
    rc  = U_k(6,:) - lambdar;
    
    % Angular acceleration rates using gradient function
    pdot = sgolayfilt(gradient(pc,dt), 3, 31);
    qdot = sgolayfilt(gradient(qc,dt), 3, 101);
    rdot = sgolayfilt(gradient(rc,dt), 3, 31);
    
    % Measured Outputs (hence still including noise)
    Vm       = Z_k(10,:);       % m denotes 'measured'
    alpham   = Z_k(11,:)';
    betam    = Z_k(12,:)';

    qnorm   = (qc./(2*Vm))'.*c;
    rnorm   = (rc./(2*Vm))'.*b;
    pnorm   = (pc./(2*Vm))'.*b;
  
    % Computed Dimensionless Aerodynamic Force Body and Moment Components 
    % Computed Using Measured Data Corrected for Biases
    OLS.ID.Y.YCx_datasets = [OLS.ID.Y.YCx_datasets; (m*Axc./(1/2*rho*Vm.^2*S))'];
    OLS.ID.Y.YCy_datasets = [OLS.ID.Y.YCy_datasets; (m*Ayc./(1/2*rho*Vm.^2*S))'];
    OLS.ID.Y.YCz_datasets = [OLS.ID.Y.YCz_datasets; (m*Azc./(1/2*rho*Vm.^2*S))'];
    OLS.ID.Y.YCl_datasets = [OLS.ID.Y.YCl_datasets; ((pdot*Ixx + qc.*rc*(Izz-Iyy) - (pc.*qc+rdot)*Ixz) ./ (1/2*rho*Vm.^2*S*b))'];
    OLS.ID.Y.YCm_datasets = [OLS.ID.Y.YCm_datasets; ((qdot*Iyy + rc.*pc*(Ixx-Izz) - (pc.^2-rc.^2)*Ixz) ./ (1/2*rho*Vm.^2*S*c))'];
    OLS.ID.Y.YCn_datasets = [OLS.ID.Y.YCn_datasets; ((rdot*Izz + pc.*qc*(Iyy-Ixx) + (qc.*rc-pdot)*Ixz) ./ (1/2*rho*Vm.^2*S*b))'];

    % OLS Parameter Estimation Original Model Structure
    % Computed Using Measured Data Corrected for Biases
    OLS.ID.X.XCx_datasets = [OLS.ID.X.XCx_datasets; [ones(N,1), alpham, alpham.^2, qnorm, de, Tc1+Tc2]];
    OLS.ID.X.XCy_datasets = [OLS.ID.X.XCy_datasets; [ones(N,1), betam, pnorm, rnorm, da, dr]];
    OLS.ID.X.XCz_datasets = [OLS.ID.X.XCz_datasets; [ones(N,1), alpham, qnorm, de, Tc1+Tc2]];
    OLS.ID.X.XCl_datasets = [OLS.ID.X.XCl_datasets; [ones(N,1), betam, pnorm, rnorm, da, dr]];
    OLS.ID.X.XCm_datasets = [OLS.ID.X.XCm_datasets; [ones(N,1), alpham, qnorm, de, Tc1+Tc2]];
    OLS.ID.X.XCn_datasets = [OLS.ID.X.XCn_datasets; [ones(N,1), betam, pnorm, rnorm, da, dr]];

    % OLS Parameter Estimation Alternative Model Structure
    % Computed Using Measured Data Corrected for Biases
    M2structure;
    OLS.ID.M2.X.XCx = [OLS.ID.M2.X.XCx; M2.Cx];
    OLS.ID.M2.X.XCy = [OLS.ID.M2.X.XCy; M2.Cy];
    OLS.ID.M2.X.XCz = [OLS.ID.M2.X.XCz; M2.Cz];
    OLS.ID.M2.X.XCl = [OLS.ID.M2.X.XCl; M2.Cl];
    OLS.ID.M2.X.XCm = [OLS.ID.M2.X.XCm; M2.Cm];
    OLS.ID.M2.X.XCn = [OLS.ID.M2.X.XCn; M2.Cn];

end

% Covariance Matrix Calculations
% Original Model Structure
OLS.ID.COV.Cx = inv((OLS.ID.X.XCx_datasets'*OLS.ID.X.XCx_datasets));
OLS.ID.COV.Cy = inv((OLS.ID.X.XCy_datasets'*OLS.ID.X.XCy_datasets));
OLS.ID.COV.Cz = inv((OLS.ID.X.XCz_datasets'*OLS.ID.X.XCz_datasets));
OLS.ID.COV.Cl = inv((OLS.ID.X.XCl_datasets'*OLS.ID.X.XCl_datasets));
OLS.ID.COV.Cm = inv((OLS.ID.X.XCm_datasets'*OLS.ID.X.XCm_datasets));
OLS.ID.COV.Cn = inv((OLS.ID.X.XCn_datasets'*OLS.ID.X.XCn_datasets));
% Alternative Model Structure
OLS.ID.M2.COV.Cx = inv((OLS.ID.M2.X.XCx'*OLS.ID.M2.X.XCx));
OLS.ID.M2.COV.Cy = inv((OLS.ID.M2.X.XCy'*OLS.ID.M2.X.XCy));
OLS.ID.M2.COV.Cz = inv((OLS.ID.M2.X.XCz'*OLS.ID.M2.X.XCz));
OLS.ID.M2.COV.Cl = inv((OLS.ID.M2.X.XCl'*OLS.ID.M2.X.XCl));
OLS.ID.M2.COV.Cm = inv((OLS.ID.M2.X.XCm'*OLS.ID.M2.X.XCm));
OLS.ID.M2.COV.Cn = inv((OLS.ID.M2.X.XCn'*OLS.ID.M2.X.XCn));

% OLS Parameter Estimation
% Original Model Structure
OLS.ID.Param.ThetaCx = OLS.ID.COV.Cx * OLS.ID.X.XCx_datasets' * OLS.ID.Y.YCx_datasets;
OLS.ID.Param.ThetaCy = OLS.ID.COV.Cy * OLS.ID.X.XCy_datasets' * OLS.ID.Y.YCy_datasets;
OLS.ID.Param.ThetaCz = OLS.ID.COV.Cz * OLS.ID.X.XCz_datasets' * OLS.ID.Y.YCz_datasets;
OLS.ID.Param.ThetaCl = OLS.ID.COV.Cl * OLS.ID.X.XCl_datasets' * OLS.ID.Y.YCl_datasets;
OLS.ID.Param.ThetaCm = OLS.ID.COV.Cm * OLS.ID.X.XCm_datasets' * OLS.ID.Y.YCm_datasets;
OLS.ID.Param.ThetaCn = OLS.ID.COV.Cn * OLS.ID.X.XCn_datasets' * OLS.ID.Y.YCn_datasets;
% Alternative Model Structure
OLS.ID.M2.Param.ThetaCx = OLS.ID.M2.COV.Cx * OLS.ID.M2.X.XCx' * OLS.ID.Y.YCx_datasets;
OLS.ID.M2.Param.ThetaCy = OLS.ID.M2.COV.Cy * OLS.ID.M2.X.XCy' * OLS.ID.Y.YCy_datasets;
OLS.ID.M2.Param.ThetaCz = OLS.ID.M2.COV.Cz * OLS.ID.M2.X.XCz' * OLS.ID.Y.YCz_datasets;
OLS.ID.M2.Param.ThetaCl = OLS.ID.M2.COV.Cl * OLS.ID.M2.X.XCl' * OLS.ID.Y.YCl_datasets;
OLS.ID.M2.Param.ThetaCm = OLS.ID.M2.COV.Cm * OLS.ID.M2.X.XCm' * OLS.ID.Y.YCm_datasets;
OLS.ID.M2.Param.ThetaCn = OLS.ID.M2.COV.Cn * OLS.ID.M2.X.XCn' * OLS.ID.Y.YCn_datasets;

% Model Calculation
% Original Model Structure
OLS.ID.Model.Cx = OLS.ID.X.XCx_datasets * OLS.ID.Param.ThetaCx;
OLS.ID.Model.Cy = OLS.ID.X.XCy_datasets * OLS.ID.Param.ThetaCy;
OLS.ID.Model.Cz = OLS.ID.X.XCz_datasets * OLS.ID.Param.ThetaCz;
OLS.ID.Model.Cl = OLS.ID.X.XCl_datasets * OLS.ID.Param.ThetaCl;
OLS.ID.Model.Cm = OLS.ID.X.XCm_datasets * OLS.ID.Param.ThetaCm;
OLS.ID.Model.Cn = OLS.ID.X.XCn_datasets * OLS.ID.Param.ThetaCn;
% Alternative Model Structure
OLS.ID.M2.Model.Cx = OLS.ID.M2.X.XCx * OLS.ID.M2.Param.ThetaCx;
OLS.ID.M2.Model.Cy = OLS.ID.M2.X.XCy * OLS.ID.M2.Param.ThetaCy;
OLS.ID.M2.Model.Cz = OLS.ID.M2.X.XCz * OLS.ID.M2.Param.ThetaCz;
OLS.ID.M2.Model.Cl = OLS.ID.M2.X.XCl * OLS.ID.M2.Param.ThetaCl;
OLS.ID.M2.Model.Cm = OLS.ID.M2.X.XCm * OLS.ID.M2.Param.ThetaCm;
OLS.ID.M2.Model.Cn = OLS.ID.M2.X.XCn * OLS.ID.M2.Param.ThetaCn;

% Residuals
% Original Model Structure
OLS.ID.Residuals.Cx = OLS.ID.Y.YCx_datasets - OLS.ID.Model.Cx;
OLS.ID.Residuals.Cy = OLS.ID.Y.YCy_datasets - OLS.ID.Model.Cy;
OLS.ID.Residuals.Cz = OLS.ID.Y.YCz_datasets - OLS.ID.Model.Cz;
OLS.ID.Residuals.Cl = OLS.ID.Y.YCl_datasets - OLS.ID.Model.Cl;
OLS.ID.Residuals.Cm = OLS.ID.Y.YCm_datasets - OLS.ID.Model.Cm;
OLS.ID.Residuals.Cn = OLS.ID.Y.YCn_datasets - OLS.ID.Model.Cn;
% Alternative Model Structure
OLS.ID.M2.Residuals.Cx = OLS.ID.Y.YCx_datasets - OLS.ID.M2.Model.Cx;
OLS.ID.M2.Residuals.Cy = OLS.ID.Y.YCy_datasets - OLS.ID.M2.Model.Cy;
OLS.ID.M2.Residuals.Cz = OLS.ID.Y.YCz_datasets - OLS.ID.M2.Model.Cz;
OLS.ID.M2.Residuals.Cl = OLS.ID.Y.YCl_datasets - OLS.ID.M2.Model.Cl;
OLS.ID.M2.Residuals.Cm = OLS.ID.Y.YCm_datasets - OLS.ID.M2.Model.Cm;
OLS.ID.M2.Residuals.Cn = OLS.ID.Y.YCn_datasets - OLS.ID.M2.Model.Cn;

% RMS Residuals
% Original Model Structure
OLS.ID.RMSResiduals.Cx = sqrt(mean(abs(OLS.ID.Residuals.Cx)).^2);
OLS.ID.RMSResiduals.Cy = sqrt(mean(abs(OLS.ID.Residuals.Cy)).^2);
OLS.ID.RMSResiduals.Cz = sqrt(mean(abs(OLS.ID.Residuals.Cz)).^2);
OLS.ID.RMSResiduals.Cl = sqrt(mean(abs(OLS.ID.Residuals.Cl)).^2);
OLS.ID.RMSResiduals.Cm = sqrt(mean(abs(OLS.ID.Residuals.Cm)).^2);
OLS.ID.RMSResiduals.Cn = sqrt(mean(abs(OLS.ID.Residuals.Cn)).^2);
% Alternative Model Structure
OLS.ID.M2.RMSResiduals.Cx = sqrt(mean(abs(OLS.ID.M2.Residuals.Cx)).^2);
OLS.ID.M2.RMSResiduals.Cy = sqrt(mean(abs(OLS.ID.M2.Residuals.Cy)).^2);
OLS.ID.M2.RMSResiduals.Cz = sqrt(mean(abs(OLS.ID.M2.Residuals.Cz)).^2);
OLS.ID.M2.RMSResiduals.Cl = sqrt(mean(abs(OLS.ID.M2.Residuals.Cl)).^2);
OLS.ID.M2.RMSResiduals.Cm = sqrt(mean(abs(OLS.ID.M2.Residuals.Cm)).^2);
OLS.ID.M2.RMSResiduals.Cn = sqrt(mean(abs(OLS.ID.M2.Residuals.Cn)).^2);

%%
% R2
% Original Model Structure
OLS.ID.R2.Cx = my_Rsquared_coeff(OLS.ID.Y.YCx_datasets, OLS.ID.Model.Cx);
OLS.ID.R2.Cy = my_Rsquared_coeff(OLS.ID.Y.YCy_datasets, OLS.ID.Model.Cy);
OLS.ID.R2.Cz = my_Rsquared_coeff(OLS.ID.Y.YCz_datasets, OLS.ID.Model.Cz);
OLS.ID.R2.Cl = my_Rsquared_coeff(OLS.ID.Y.YCl_datasets, OLS.ID.Model.Cl);
OLS.ID.R2.Cm = my_Rsquared_coeff(OLS.ID.Y.YCm_datasets, OLS.ID.Model.Cm);
OLS.ID.R2.Cn = my_Rsquared_coeff(OLS.ID.Y.YCn_datasets, OLS.ID.Model.Cn);
% Alternative Model Structure
OLS.ID.M2.R2.Cx = my_Rsquared_coeff(OLS.ID.Y.YCx_datasets, OLS.ID.M2.Model.Cx);
OLS.ID.M2.R2.Cy = my_Rsquared_coeff(OLS.ID.Y.YCy_datasets, OLS.ID.M2.Model.Cy);
OLS.ID.M2.R2.Cz = my_Rsquared_coeff(OLS.ID.Y.YCz_datasets, OLS.ID.M2.Model.Cz);
OLS.ID.M2.R2.Cl = my_Rsquared_coeff(OLS.ID.Y.YCl_datasets, OLS.ID.M2.Model.Cl);
OLS.ID.M2.R2.Cm = my_Rsquared_coeff(OLS.ID.Y.YCm_datasets, OLS.ID.M2.Model.Cm);
OLS.ID.M2.R2.Cn = my_Rsquared_coeff(OLS.ID.Y.YCn_datasets, OLS.ID.M2.Model.Cn);

disp("OLS Identification Complete")




