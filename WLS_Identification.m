%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   WLS Identification Model
%
%   Author: Y.J.E. Prencipe, based on C.C. de Visser
%   Student Number: 4777158
%   Course: AE4320 System Identification of Aerospace Vehicles
%   Place: Delft University of Technology, 2023
%   Email: y.j.e.prencipe@student.tudelft.nl
%   Version: 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate Weighing Matrix For Identification Datasets
w1 = 1./OLS.ID.Residuals.Cx.^2;
w2 = 1./OLS.ID.Residuals.Cy.^2;
w3 = 1./OLS.ID.Residuals.Cz.^2;
w4 = 1./OLS.ID.Residuals.Cl.^2;
w5 = 1./OLS.ID.Residuals.Cm.^2;
w6 = 1./OLS.ID.Residuals.Cn.^2;
% alternative model weights
w1a = 1./OLS.ID.M2.Residuals.Cx.^2;
w2a = 1./OLS.ID.M2.Residuals.Cy.^2;
w3a = 1./OLS.ID.M2.Residuals.Cz.^2;
w4a = 1./OLS.ID.M2.Residuals.Cl.^2;
w5a = 1./OLS.ID.M2.Residuals.Cm.^2;
w6a = 1./OLS.ID.M2.Residuals.Cn.^2;

N = size(w1,1);

WLS.ID.W.WCx = w1 .* speye(N,N);
WLS.ID.W.WCy = w2 .* speye(N,N);
WLS.ID.W.WCz = w3 .* speye(N,N);
WLS.ID.W.WCl = w4 .* speye(N,N);
WLS.ID.W.WCm = w5 .* speye(N,N);
WLS.ID.W.WCn = w6 .* speye(N,N);
% alternative model weighting matrices
WLS.ID.M2.W.WCx = w1a .* speye(N,N);
WLS.ID.M2.W.WCy = w2a .* speye(N,N);
WLS.ID.M2.W.WCz = w3a .* speye(N,N);
WLS.ID.M2.W.WCl = w4a .* speye(N,N);
WLS.ID.M2.W.WCm = w5a .* speye(N,N);
WLS.ID.M2.W.WCn = w6a .* speye(N,N);


% Define Concatenated X Matrices
WLS.ID.X.XCx_datasets = [];
WLS.ID.X.XCy_datasets = [];
WLS.ID.X.XCz_datasets = [];
WLS.ID.X.XCl_datasets = [];
WLS.ID.X.XCm_datasets = [];
WLS.ID.X.XCn_datasets = [];
% Define Concatenated X Matrices for Alternative Model Structure
WLS.ID.M2.X.XCx = [];
WLS.ID.M2.X.XCy = [];
WLS.ID.M2.X.XCz = [];
WLS.ID.M2.X.XCl = [];
WLS.ID.M2.X.XCm = [];
WLS.ID.M2.X.XCn = [];
% Define Concatenated Y Matrices
WLS.ID.Y.YCx_datasets = [];
WLS.ID.Y.YCy_datasets = [];
WLS.ID.Y.YCz_datasets = [];
WLS.ID.Y.YCl_datasets = [];
WLS.ID.Y.YCm_datasets = [];
WLS.ID.Y.YCn_datasets = [];

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

    qnorm     = (qc./Vm)'.*c;
    rnorm   = (rc./(2*Vm))'.*b;
    pnorm   = (pc./(2*Vm))'.*b;

    WLS.ID.Y.YCx_datasets = [WLS.ID.Y.YCx_datasets; (m*Axc./(1/2*rho*Vm.^2*S))'];
    WLS.ID.Y.YCy_datasets = [WLS.ID.Y.YCy_datasets; (m*Ayc./(1/2*rho*Vm.^2*S))'];
    WLS.ID.Y.YCz_datasets = [WLS.ID.Y.YCz_datasets; (m*Azc./(1/2*rho*Vm.^2*S))'];
    WLS.ID.Y.YCl_datasets = [WLS.ID.Y.YCl_datasets; ((pdot*Ixx + qc.*rc*(Izz-Iyy) - (pc.*qc+rdot)*Ixz) ./ (1/2*rho*Vm.^2*S*b))'];
    WLS.ID.Y.YCm_datasets = [WLS.ID.Y.YCm_datasets; ((qdot*Iyy + rc.*pc*(Ixx-Izz) - (pc.^2-rc.^2)*Ixz) ./ (1/2*rho*Vm.^2*S*c))'];
    WLS.ID.Y.YCn_datasets = [WLS.ID.Y.YCn_datasets; ((rdot*Izz + pc.*qc*(Iyy-Ixx) + (qc.*rc-pdot)*Ixz) ./ (1/2*rho*Vm.^2*S*b))'];
    
    WLS.ID.X.XCx_datasets = [WLS.ID.X.XCx_datasets; [ones(N,1), alpham, alpham.^2, qnorm, de, Tc1+Tc2]];
    WLS.ID.X.XCy_datasets = [WLS.ID.X.XCy_datasets; [ones(N,1), betam, (pc./(2*Vm))'.*b, (rc./(2*Vm))'.*b, da, dr]];
    WLS.ID.X.XCz_datasets = [WLS.ID.X.XCz_datasets; [ones(N,1), alpham, qnorm, de, Tc1+Tc2]];
    WLS.ID.X.XCl_datasets = [WLS.ID.X.XCl_datasets; [ones(N,1), betam, (pc./(2*Vm))'.*b, (rc./(2*Vm))'.*b, da, dr]];
    WLS.ID.X.XCm_datasets = [WLS.ID.X.XCm_datasets; [ones(N,1), alpham, qnorm, de, Tc1+Tc2]];
    WLS.ID.X.XCn_datasets = [WLS.ID.X.XCn_datasets; [ones(N,1), betam, (pc./(2*Vm))'.*b, (rc./(2*Vm))'.*b, da, dr]];

    % WLS Parameter Estimation Alternative Model Structure
    % Computed Using Measured Data Corrected for Biases
    M2structure;
    WLS.ID.M2.X.XCx = [WLS.ID.M2.X.XCx; M2.Cx];
    WLS.ID.M2.X.XCy = [WLS.ID.M2.X.XCy; M2.Cy];
    WLS.ID.M2.X.XCz = [WLS.ID.M2.X.XCz; M2.Cz];
    WLS.ID.M2.X.XCl = [WLS.ID.M2.X.XCl; M2.Cl];
    WLS.ID.M2.X.XCm = [WLS.ID.M2.X.XCm; M2.Cm];
    WLS.ID.M2.X.XCn = [WLS.ID.M2.X.XCn; M2.Cn];

end

% Covariance Matrix Calculations
WLS.ID.COV.Cx = inv(WLS.ID.X.XCx_datasets' * (WLS.ID.W.WCx) * WLS.ID.X.XCx_datasets);
WLS.ID.COV.Cy = inv(WLS.ID.X.XCy_datasets' * (WLS.ID.W.WCy) * WLS.ID.X.XCy_datasets);
WLS.ID.COV.Cz = inv(WLS.ID.X.XCz_datasets' * (WLS.ID.W.WCz) * WLS.ID.X.XCz_datasets);
WLS.ID.COV.Cl = inv(WLS.ID.X.XCl_datasets' * (WLS.ID.W.WCl) * WLS.ID.X.XCl_datasets);
WLS.ID.COV.Cm = inv(WLS.ID.X.XCm_datasets' * (WLS.ID.W.WCm) * WLS.ID.X.XCm_datasets);
WLS.ID.COV.Cn = inv(WLS.ID.X.XCn_datasets' * (WLS.ID.W.WCn) * WLS.ID.X.XCn_datasets);
% Alternative Model Structure
WLS.ID.M2.COV.Cx = inv((WLS.ID.M2.X.XCx'* (WLS.ID.M2.W.WCx) * WLS.ID.M2.X.XCx));
WLS.ID.M2.COV.Cy = inv((WLS.ID.M2.X.XCy'* (WLS.ID.M2.W.WCy) * WLS.ID.M2.X.XCy));
WLS.ID.M2.COV.Cz = inv((WLS.ID.M2.X.XCz'* (WLS.ID.M2.W.WCz) * WLS.ID.M2.X.XCz));
WLS.ID.M2.COV.Cl = inv((WLS.ID.M2.X.XCl'* (WLS.ID.M2.W.WCl) * WLS.ID.M2.X.XCl));
WLS.ID.M2.COV.Cm = inv((WLS.ID.M2.X.XCm'* (WLS.ID.M2.W.WCm) * WLS.ID.M2.X.XCm));
WLS.ID.M2.COV.Cn = inv((WLS.ID.M2.X.XCn'* (WLS.ID.M2.W.WCn) * WLS.ID.M2.X.XCn));

% WLS Parameter Estimation
WLS.ID.Param.ThetaCx =  WLS.ID.COV.Cx * WLS.ID.X.XCx_datasets' * (WLS.ID.W.WCx) * WLS.ID.Y.YCx_datasets;
WLS.ID.Param.ThetaCy =  WLS.ID.COV.Cy * WLS.ID.X.XCy_datasets' * (WLS.ID.W.WCy) * WLS.ID.Y.YCy_datasets;
WLS.ID.Param.ThetaCz =  WLS.ID.COV.Cz * WLS.ID.X.XCz_datasets' * (WLS.ID.W.WCz) * WLS.ID.Y.YCz_datasets;
WLS.ID.Param.ThetaCl =  WLS.ID.COV.Cl * WLS.ID.X.XCl_datasets' * (WLS.ID.W.WCl) * WLS.ID.Y.YCl_datasets;
WLS.ID.Param.ThetaCm =  WLS.ID.COV.Cm * WLS.ID.X.XCm_datasets' * (WLS.ID.W.WCm) * WLS.ID.Y.YCm_datasets;
WLS.ID.Param.ThetaCn =  WLS.ID.COV.Cn * WLS.ID.X.XCn_datasets' * (WLS.ID.W.WCn) * WLS.ID.Y.YCn_datasets;
% Alternative Model Structure
WLS.ID.M2.Param.ThetaCx = WLS.ID.M2.COV.Cx * WLS.ID.M2.X.XCx' * (WLS.ID.M2.W.WCx) * WLS.ID.Y.YCx_datasets;
WLS.ID.M2.Param.ThetaCy = WLS.ID.M2.COV.Cy * WLS.ID.M2.X.XCy' * (WLS.ID.M2.W.WCy) * WLS.ID.Y.YCy_datasets;
WLS.ID.M2.Param.ThetaCz = WLS.ID.M2.COV.Cz * WLS.ID.M2.X.XCz' * (WLS.ID.M2.W.WCz) * WLS.ID.Y.YCz_datasets;
WLS.ID.M2.Param.ThetaCl = WLS.ID.M2.COV.Cl * WLS.ID.M2.X.XCl' * (WLS.ID.M2.W.WCl) * WLS.ID.Y.YCl_datasets;
WLS.ID.M2.Param.ThetaCm = WLS.ID.M2.COV.Cm * WLS.ID.M2.X.XCm' * (WLS.ID.M2.W.WCm) * WLS.ID.Y.YCm_datasets;
WLS.ID.M2.Param.ThetaCn = WLS.ID.M2.COV.Cn * WLS.ID.M2.X.XCn' * (WLS.ID.M2.W.WCn) * WLS.ID.Y.YCn_datasets;

% Model Calculation
WLS.ID.Model.Cx = WLS.ID.X.XCx_datasets * WLS.ID.Param.ThetaCx;
WLS.ID.Model.Cy = WLS.ID.X.XCy_datasets * WLS.ID.Param.ThetaCy;
WLS.ID.Model.Cz = WLS.ID.X.XCz_datasets * WLS.ID.Param.ThetaCz;
WLS.ID.Model.Cl = WLS.ID.X.XCl_datasets * WLS.ID.Param.ThetaCl;
WLS.ID.Model.Cm = WLS.ID.X.XCm_datasets * WLS.ID.Param.ThetaCm;
WLS.ID.Model.Cn = WLS.ID.X.XCn_datasets * WLS.ID.Param.ThetaCn;
% Alternative Model Structure
WLS.ID.M2.Model.Cx = WLS.ID.M2.X.XCx * WLS.ID.M2.Param.ThetaCx;
WLS.ID.M2.Model.Cy = WLS.ID.M2.X.XCy * WLS.ID.M2.Param.ThetaCy;
WLS.ID.M2.Model.Cz = WLS.ID.M2.X.XCz * WLS.ID.M2.Param.ThetaCz;
WLS.ID.M2.Model.Cl = WLS.ID.M2.X.XCl * WLS.ID.M2.Param.ThetaCl;
WLS.ID.M2.Model.Cm = WLS.ID.M2.X.XCm * WLS.ID.M2.Param.ThetaCm;
WLS.ID.M2.Model.Cn = WLS.ID.M2.X.XCn * WLS.ID.M2.Param.ThetaCn;

% Residuals
WLS.ID.Residuals.Cx = WLS.ID.Y.YCx_datasets - WLS.ID.Model.Cx;
WLS.ID.Residuals.Cy = WLS.ID.Y.YCy_datasets - WLS.ID.Model.Cy;
WLS.ID.Residuals.Cz = WLS.ID.Y.YCz_datasets - WLS.ID.Model.Cz;
WLS.ID.Residuals.Cl = WLS.ID.Y.YCl_datasets - WLS.ID.Model.Cl;
WLS.ID.Residuals.Cm = WLS.ID.Y.YCm_datasets - WLS.ID.Model.Cm;
WLS.ID.Residuals.Cn = WLS.ID.Y.YCn_datasets - WLS.ID.Model.Cn;
% Alternative Model Structure
WLS.ID.M2.Residuals.Cx = WLS.ID.Y.YCx_datasets - WLS.ID.M2.Model.Cx;
WLS.ID.M2.Residuals.Cy = WLS.ID.Y.YCy_datasets - WLS.ID.M2.Model.Cy;
WLS.ID.M2.Residuals.Cz = WLS.ID.Y.YCz_datasets - WLS.ID.M2.Model.Cz;
WLS.ID.M2.Residuals.Cl = WLS.ID.Y.YCl_datasets - WLS.ID.M2.Model.Cl;
WLS.ID.M2.Residuals.Cm = WLS.ID.Y.YCm_datasets - WLS.ID.M2.Model.Cm;
WLS.ID.M2.Residuals.Cn = WLS.ID.Y.YCn_datasets - WLS.ID.M2.Model.Cn;

% RMS Residuals
WLS.ID.RMSResiduals.Cx = sqrt(mean(abs(WLS.ID.Residuals.Cx)).^2);
WLS.ID.RMSResiduals.Cy = sqrt(mean(abs(WLS.ID.Residuals.Cy)).^2);
WLS.ID.RMSResiduals.Cz = sqrt(mean(abs(WLS.ID.Residuals.Cz)).^2);
WLS.ID.RMSResiduals.Cl = sqrt(mean(abs(WLS.ID.Residuals.Cl)).^2);
WLS.ID.RMSResiduals.Cm = sqrt(mean(abs(WLS.ID.Residuals.Cm)).^2);
WLS.ID.RMSResiduals.Cn = sqrt(mean(abs(WLS.ID.Residuals.Cn)).^2);
% Alternative Model Structure
WLS.ID.M2.RMSResiduals.Cx = sqrt(mean(abs(WLS.ID.M2.Residuals.Cx)).^2);
WLS.ID.M2.RMSResiduals.Cy = sqrt(mean(abs(WLS.ID.M2.Residuals.Cy)).^2);
WLS.ID.M2.RMSResiduals.Cz = sqrt(mean(abs(WLS.ID.M2.Residuals.Cz)).^2);
WLS.ID.M2.RMSResiduals.Cl = sqrt(mean(abs(WLS.ID.M2.Residuals.Cl)).^2);
WLS.ID.M2.RMSResiduals.Cm = sqrt(mean(abs(WLS.ID.M2.Residuals.Cm)).^2);
WLS.ID.M2.RMSResiduals.Cn = sqrt(mean(abs(WLS.ID.M2.Residuals.Cn)).^2);

%%
% R2
WLS.ID.R2.Cx = my_Rsquared_coeff(WLS.ID.Y.YCx_datasets, WLS.ID.Model.Cx);
WLS.ID.R2.Cy = my_Rsquared_coeff(WLS.ID.Y.YCy_datasets, WLS.ID.Model.Cy);
WLS.ID.R2.Cz = my_Rsquared_coeff(WLS.ID.Y.YCz_datasets, WLS.ID.Model.Cz);
WLS.ID.R2.Cl = my_Rsquared_coeff(WLS.ID.Y.YCl_datasets, WLS.ID.Model.Cl);
WLS.ID.R2.Cm = my_Rsquared_coeff(WLS.ID.Y.YCm_datasets, WLS.ID.Model.Cm);
WLS.ID.R2.Cn = my_Rsquared_coeff(WLS.ID.Y.YCn_datasets, WLS.ID.Model.Cn);
% Alternative Model Structure
WLS.ID.M2.R2.Cx = my_Rsquared_coeff(WLS.ID.Y.YCx_datasets, WLS.ID.M2.Model.Cx);
WLS.ID.M2.R2.Cy = my_Rsquared_coeff(WLS.ID.Y.YCy_datasets, WLS.ID.M2.Model.Cy);
WLS.ID.M2.R2.Cz = my_Rsquared_coeff(WLS.ID.Y.YCz_datasets, WLS.ID.M2.Model.Cz);
WLS.ID.M2.R2.Cl = my_Rsquared_coeff(WLS.ID.Y.YCl_datasets, WLS.ID.M2.Model.Cl);
WLS.ID.M2.R2.Cm = my_Rsquared_coeff(WLS.ID.Y.YCm_datasets, WLS.ID.M2.Model.Cm);
WLS.ID.M2.R2.Cn = my_Rsquared_coeff(WLS.ID.Y.YCn_datasets, WLS.ID.M2.Model.Cn);


disp("WLS Identification Complete")





