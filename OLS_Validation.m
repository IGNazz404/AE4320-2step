%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OLS Validation Model
%
%   Author: Y.J.E. Prencipe, based on C.C. de Visser
%   Student Number: 4777158
%   Course: AE4320 System Identification of Aerospace Vehicles
%   Place: Delft University of Technology, 2023
%   Email: y.j.e.prencipe@student.tudelft.nl
%   Version: 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Concatenated X Matrices
OLS.VAL.X.XCx_datasets = [];
OLS.VAL.X.XCy_datasets = [];
OLS.VAL.X.XCz_datasets = [];
OLS.VAL.X.XCl_datasets = [];
OLS.VAL.X.XCm_datasets = [];
OLS.VAL.X.XCn_datasets = [];
% Define Concatenated X Matrices for Alternative Model Structure
OLS.VAL.M2.X.XCx = [];
OLS.VAL.M2.X.XCy = [];
OLS.VAL.M2.X.XCz = [];
OLS.VAL.M2.X.XCl = [];
OLS.VAL.M2.X.XCm = [];
OLS.VAL.M2.X.XCn = [];
% Define Concatenated Y Matrices
OLS.VAL.Y.YCx_datasets = [];
OLS.VAL.Y.YCy_datasets = [];
OLS.VAL.Y.YCz_datasets = [];
OLS.VAL.Y.YCl_datasets = [];
OLS.VAL.Y.YCm_datasets = [];
OLS.VAL.Y.YCn_datasets = [];

for i=5:6   % CHOSE THESE 2 BECAUSE THEY HAVE FULL MODEL COVERAGE - I.E. EVERY MODEL PARAMETER HAS SOME KIND OF EXCITEMENT USING THESE MANOEUVRES

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
    OLS.VAL.Y.YCx_datasets = [OLS.VAL.Y.YCx_datasets; (m*Axc./(1/2*rho*Vm.^2*S))'];
    OLS.VAL.Y.YCy_datasets = [OLS.VAL.Y.YCy_datasets; (m*Ayc./(1/2*rho*Vm.^2*S))'];
    OLS.VAL.Y.YCz_datasets = [OLS.VAL.Y.YCz_datasets; (m*Azc./(1/2*rho*Vm.^2*S))'];
    OLS.VAL.Y.YCl_datasets = [OLS.VAL.Y.YCl_datasets; ((pdot*Ixx + qc.*rc*(Izz-Iyy) - (pc.*qc+rdot)*Ixz) ./ (1/2*rho*Vm.^2*S*b))'];
    OLS.VAL.Y.YCm_datasets = [OLS.VAL.Y.YCm_datasets; ((qdot*Iyy + rc.*pc*(Ixx-Izz) - (pc.^2-rc.^2)*Ixz) ./ (1/2*rho*Vm.^2*S*c))'];
    OLS.VAL.Y.YCn_datasets = [OLS.VAL.Y.YCn_datasets; ((rdot*Izz + pc.*qc*(Iyy-Ixx) + (qc.*rc-pdot)*Ixz) ./ (1/2*rho*Vm.^2*S*b))'];

    % OLS Parameter Estimation
    % Computed Using Measured Data Corrected for Biases
    OLS.VAL.X.XCx_datasets = [OLS.VAL.X.XCx_datasets; [ones(N,1), alpham, alpham.^2, qnorm, de, Tc1+Tc2]];
    OLS.VAL.X.XCy_datasets = [OLS.VAL.X.XCy_datasets; [ones(N,1), betam, (pc./(2*Vm))'.*b, (rc./(2*Vm))'.*b, da, dr]];
    OLS.VAL.X.XCz_datasets = [OLS.VAL.X.XCz_datasets; [ones(N,1), alpham, qnorm, de, Tc1+Tc2]];
    OLS.VAL.X.XCl_datasets = [OLS.VAL.X.XCl_datasets; [ones(N,1), betam, (pc./(2*Vm))'.*b, (rc./(2*Vm))'.*b, da, dr]];
    OLS.VAL.X.XCm_datasets = [OLS.VAL.X.XCm_datasets; [ones(N,1), alpham, qnorm, de, Tc1+Tc2]];
    OLS.VAL.X.XCn_datasets = [OLS.VAL.X.XCn_datasets; [ones(N,1), betam, (pc./(2*Vm))'.*b, (rc./(2*Vm))'.*b, da, dr]];

    M2structure;
    OLS.VAL.M2.X.XCx = [OLS.VAL.M2.X.XCx; M2.Cx];
    OLS.VAL.M2.X.XCy = [OLS.VAL.M2.X.XCy; M2.Cy];
    OLS.VAL.M2.X.XCz = [OLS.VAL.M2.X.XCz; M2.Cz];
    OLS.VAL.M2.X.XCl = [OLS.VAL.M2.X.XCl; M2.Cl];
    OLS.VAL.M2.X.XCm = [OLS.VAL.M2.X.XCm; M2.Cm];
    OLS.VAL.M2.X.XCn = [OLS.VAL.M2.X.XCn; M2.Cn];

end

% Model Calculation
OLS.VAL.Model.Cx = OLS.VAL.X.XCx_datasets * OLS.ID.Param.ThetaCx;     % !!!! USES IDENTIFIED PARAMETERS !!!! %
OLS.VAL.Model.Cy = OLS.VAL.X.XCy_datasets * OLS.ID.Param.ThetaCy;
OLS.VAL.Model.Cz = OLS.VAL.X.XCz_datasets * OLS.ID.Param.ThetaCz;
OLS.VAL.Model.Cl = OLS.VAL.X.XCl_datasets * OLS.ID.Param.ThetaCl;
OLS.VAL.Model.Cm = OLS.VAL.X.XCm_datasets * OLS.ID.Param.ThetaCm;
OLS.VAL.Model.Cn = OLS.VAL.X.XCn_datasets * OLS.ID.Param.ThetaCn;
% Alternative Model Structure
OLS.VAL.M2.Model.Cx = OLS.VAL.M2.X.XCx * OLS.ID.M2.Param.ThetaCx;
OLS.VAL.M2.Model.Cy = OLS.VAL.M2.X.XCy * OLS.ID.M2.Param.ThetaCy;
OLS.VAL.M2.Model.Cz = OLS.VAL.M2.X.XCz * OLS.ID.M2.Param.ThetaCz;
OLS.VAL.M2.Model.Cl = OLS.VAL.M2.X.XCl * OLS.ID.M2.Param.ThetaCl;
OLS.VAL.M2.Model.Cm = OLS.VAL.M2.X.XCm * OLS.ID.M2.Param.ThetaCm;
OLS.VAL.M2.Model.Cn = OLS.VAL.M2.X.XCn * OLS.ID.M2.Param.ThetaCn;

% Residuals
OLS.VAL.Residuals.Cx = OLS.VAL.Y.YCx_datasets - OLS.VAL.Model.Cx;
OLS.VAL.Residuals.Cy = OLS.VAL.Y.YCy_datasets - OLS.VAL.Model.Cy;
OLS.VAL.Residuals.Cz = OLS.VAL.Y.YCz_datasets - OLS.VAL.Model.Cz;
OLS.VAL.Residuals.Cl = OLS.VAL.Y.YCl_datasets - OLS.VAL.Model.Cl;
OLS.VAL.Residuals.Cm = OLS.VAL.Y.YCm_datasets - OLS.VAL.Model.Cm;
OLS.VAL.Residuals.Cn = OLS.VAL.Y.YCn_datasets - OLS.VAL.Model.Cn;
% Alternative Model Structure
OLS.VAL.M2.Residuals.Cx = OLS.VAL.Y.YCx_datasets - OLS.VAL.M2.Model.Cx;
OLS.VAL.M2.Residuals.Cy = OLS.VAL.Y.YCy_datasets - OLS.VAL.M2.Model.Cy;
OLS.VAL.M2.Residuals.Cz = OLS.VAL.Y.YCz_datasets - OLS.VAL.M2.Model.Cz;
OLS.VAL.M2.Residuals.Cl = OLS.VAL.Y.YCl_datasets - OLS.VAL.M2.Model.Cl;
OLS.VAL.M2.Residuals.Cm = OLS.VAL.Y.YCm_datasets - OLS.VAL.M2.Model.Cm;
OLS.VAL.M2.Residuals.Cn = OLS.VAL.Y.YCn_datasets - OLS.VAL.M2.Model.Cn;

% RMS Residuals
OLS.VAL.RMSResiduals.Cx = sqrt(mean(abs(OLS.VAL.Residuals.Cx)).^2);
OLS.VAL.RMSResiduals.Cy = sqrt(mean(abs(OLS.VAL.Residuals.Cy)).^2);
OLS.VAL.RMSResiduals.Cz = sqrt(mean(abs(OLS.VAL.Residuals.Cz)).^2);
OLS.VAL.RMSResiduals.Cl = sqrt(mean(abs(OLS.VAL.Residuals.Cl)).^2);
OLS.VAL.RMSResiduals.Cm = sqrt(mean(abs(OLS.VAL.Residuals.Cm)).^2);
OLS.VAL.RMSResiduals.Cn = sqrt(mean(abs(OLS.VAL.Residuals.Cn)).^2);
% Alternative Model Structure
OLS.VAL.M2.RMSResiduals.Cx = sqrt(mean(abs(OLS.VAL.M2.Residuals.Cx)).^2);
OLS.VAL.M2.RMSResiduals.Cy = sqrt(mean(abs(OLS.VAL.M2.Residuals.Cy)).^2);
OLS.VAL.M2.RMSResiduals.Cz = sqrt(mean(abs(OLS.VAL.M2.Residuals.Cz)).^2);
OLS.VAL.M2.RMSResiduals.Cl = sqrt(mean(abs(OLS.VAL.M2.Residuals.Cl)).^2);
OLS.VAL.M2.RMSResiduals.Cm = sqrt(mean(abs(OLS.VAL.M2.Residuals.Cm)).^2);
OLS.VAL.M2.RMSResiduals.Cn = sqrt(mean(abs(OLS.VAL.M2.Residuals.Cn)).^2);

%%
% R2
OLS.VAL.R2.Cx = my_Rsquared_coeff(OLS.VAL.Y.YCx_datasets, OLS.VAL.Model.Cx);
OLS.VAL.R2.Cy = my_Rsquared_coeff(OLS.VAL.Y.YCy_datasets, OLS.VAL.Model.Cy);
OLS.VAL.R2.Cz = my_Rsquared_coeff(OLS.VAL.Y.YCz_datasets, OLS.VAL.Model.Cz);
OLS.VAL.R2.Cl = my_Rsquared_coeff(OLS.VAL.Y.YCl_datasets, OLS.VAL.Model.Cl);
OLS.VAL.R2.Cm = my_Rsquared_coeff(OLS.VAL.Y.YCm_datasets, OLS.VAL.Model.Cm);
OLS.VAL.R2.Cn = my_Rsquared_coeff(OLS.VAL.Y.YCn_datasets, OLS.VAL.Model.Cn);
% Alternative Model Structure
OLS.VAL.M2.R2.Cx = my_Rsquared_coeff(OLS.VAL.Y.YCx_datasets, OLS.VAL.M2.Model.Cx);
OLS.VAL.M2.R2.Cy = my_Rsquared_coeff(OLS.VAL.Y.YCy_datasets, OLS.VAL.M2.Model.Cy);
OLS.VAL.M2.R2.Cz = my_Rsquared_coeff(OLS.VAL.Y.YCz_datasets, OLS.VAL.M2.Model.Cz);
OLS.VAL.M2.R2.Cl = my_Rsquared_coeff(OLS.VAL.Y.YCl_datasets, OLS.VAL.M2.Model.Cl);
OLS.VAL.M2.R2.Cm = my_Rsquared_coeff(OLS.VAL.Y.YCm_datasets, OLS.VAL.M2.Model.Cm);
OLS.VAL.M2.R2.Cn = my_Rsquared_coeff(OLS.VAL.Y.YCn_datasets, OLS.VAL.M2.Model.Cn);


disp("OLS Validation Complete")
