%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Model Term Dominance Tester for Least Square Model Analysis
%
%   Author: Y.J.E. Prencipe
%   Student Number: 4777158
%   Course: AE4320 System Identification of Aerospace Vehicles
%   Place: Delft University of Technology, 2023
%   Email: y.j.e.prencipe@student.tudelft.nl
%   Version: 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First run paramaterEstimation.m

numDataPoints = 100;
dom.numVariables.Cx = length(OLS.ID.Param.ThetaCx);
dom.numVariables.Cy = length(OLS.ID.Param.ThetaCy);
dom.numVariables.Cz = length(OLS.ID.Param.ThetaCz);
dom.numVariables.Cl = length(OLS.ID.Param.ThetaCl);
dom.numVariables.Cm = length(OLS.ID.Param.ThetaCm);
dom.numVariables.Cn = length(OLS.ID.Param.ThetaCn);

dom.meansOLS.Cx = OLS.ID.Param.ThetaCx(1:dom.numVariables.Cx, 1)'; 
dom.meansOLS.Cy = OLS.ID.Param.ThetaCy(1:dom.numVariables.Cy, 1)'; 
dom.meansOLS.Cz = OLS.ID.Param.ThetaCz(1:dom.numVariables.Cz, 1)'; 
dom.meansOLS.Cl = OLS.ID.Param.ThetaCl(1:dom.numVariables.Cl, 1)'; 
dom.meansOLS.Cm = OLS.ID.Param.ThetaCm(1:dom.numVariables.Cm, 1)'; 
dom.meansOLS.Cn = OLS.ID.Param.ThetaCn(1:dom.numVariables.Cn, 1)'; 

dom.variancesOLS.Cx = 10.^floor(log10(abs(OLS.ID.Param.ThetaCx(1:dom.numVariables.Cx, 1)))); dom.variancesOLS.Cx = dom.variancesOLS.Cx';
dom.variancesOLS.Cy = 10.^floor(log10(abs(OLS.ID.Param.ThetaCy(1:dom.numVariables.Cy, 1)))); dom.variancesOLS.Cy = dom.variancesOLS.Cy';
dom.variancesOLS.Cz = 10.^floor(log10(abs(OLS.ID.Param.ThetaCz(1:dom.numVariables.Cz, 1)))); dom.variancesOLS.Cz = dom.variancesOLS.Cz';
dom.variancesOLS.Cl = 10.^floor(log10(abs(OLS.ID.Param.ThetaCl(1:dom.numVariables.Cl, 1)))); dom.variancesOLS.Cl = dom.variancesOLS.Cl';
dom.variancesOLS.Cm = 10.^floor(log10(abs(OLS.ID.Param.ThetaCm(1:dom.numVariables.Cm, 1)))); dom.variancesOLS.Cm = dom.variancesOLS.Cm';
dom.variancesOLS.Cn = 10.^floor(log10(abs(OLS.ID.Param.ThetaCn(1:dom.numVariables.Cn, 1)))); dom.variancesOLS.Cn = dom.variancesOLS.Cn';

dom.dataMatrixOLS.Cx = zeros(numDataPoints, dom.numVariables.Cx);
dom.dataMatrixOLS.Cy = zeros(numDataPoints, dom.numVariables.Cy);
dom.dataMatrixOLS.Cz = zeros(numDataPoints, dom.numVariables.Cz);
dom.dataMatrixOLS.Cl = zeros(numDataPoints, dom.numVariables.Cl);
dom.dataMatrixOLS.Cm = zeros(numDataPoints, dom.numVariables.Cm);
dom.dataMatrixOLS.Cn = zeros(numDataPoints, dom.numVariables.Cn);

dom.rmsMatrixOLS.Cx = zeros(numDataPoints, dom.numVariables.Cx);
dom.rmsMatrixOLS.Cy = zeros(numDataPoints, dom.numVariables.Cy);
dom.rmsMatrixOLS.Cz = zeros(numDataPoints, dom.numVariables.Cz);
dom.rmsMatrixOLS.Cl = zeros(numDataPoints, dom.numVariables.Cl);
dom.rmsMatrixOLS.Cm = zeros(numDataPoints, dom.numVariables.Cm);
dom.rmsMatrixOLS.Cn = zeros(numDataPoints, dom.numVariables.Cn);

for i = 1:dom.numVariables.Cx
    dom.dataMatrixOLS.Cx(:, i) = normrnd(dom.meansOLS.Cx(i), sqrt(dom.variancesOLS.Cx(i)), numDataPoints, 1);
    dom.dataMatrixOLS.Cy(:, i) = normrnd(dom.meansOLS.Cy(i), sqrt(dom.variancesOLS.Cy(i)), numDataPoints, 1);
    dom.dataMatrixOLS.Cl(:, i) = normrnd(dom.meansOLS.Cl(i), sqrt(dom.variancesOLS.Cl(i)), numDataPoints, 1);
    dom.dataMatrixOLS.Cn(:, i) = normrnd(dom.meansOLS.Cn(i), sqrt(dom.variancesOLS.Cn(i)), numDataPoints, 1);
end

for i = 1:dom.numVariables.Cz
    dom.dataMatrixOLS.Cz(:, i) = normrnd(dom.meansOLS.Cz(i), sqrt(dom.variancesOLS.Cz(i)), numDataPoints, 1);
    dom.dataMatrixOLS.Cm(:, i) = normrnd(dom.meansOLS.Cm(i), sqrt(dom.variancesOLS.Cm(i)), numDataPoints, 1);
end

for i = 1:dom.numVariables.Cx
    dom.ParamOLS.Cx = OLS.ID.Param.ThetaCx;
    dom.ParamOLS.Cy = OLS.ID.Param.ThetaCy;
    dom.ParamOLS.Cl = OLS.ID.Param.ThetaCl;
    dom.ParamOLS.Cn = OLS.ID.Param.ThetaCn;

    for j = 1:numDataPoints
        % OLS
        dom.ParamOLS.Cx(i) = dom.dataMatrixOLS.Cx(j, i);
        dom.ParamOLS.Cy(i) = dom.dataMatrixOLS.Cy(j, i);
        dom.ParamOLS.Cl(i) = dom.dataMatrixOLS.Cl(j, i);
        dom.ParamOLS.Cn(i) = dom.dataMatrixOLS.Cn(j, i);

        dom.domModel.Cx = OLS.ID.X.XCx_datasets * dom.ParamOLS.Cx;
        dom.domModel.Cy = OLS.ID.X.XCy_datasets * dom.ParamOLS.Cy;
        dom.domModel.Cl = OLS.ID.X.XCl_datasets * dom.ParamOLS.Cl;
        dom.domModel.Cn = OLS.ID.X.XCn_datasets * dom.ParamOLS.Cn;

        dom.domResid.Cx = OLS.ID.Y.YCx_datasets - dom.domModel.Cx;
        dom.domResid.Cy = OLS.ID.Y.YCy_datasets - dom.domModel.Cy;
        dom.domResid.Cl = OLS.ID.Y.YCl_datasets - dom.domModel.Cl;
        dom.domResid.Cn = OLS.ID.Y.YCn_datasets - dom.domModel.Cn;

        dom.domRMSRes.Cx = sqrt(mean(abs(dom.domResid.Cx)).^2);
        dom.domRMSRes.Cy = sqrt(mean(abs(dom.domResid.Cy)).^2);
        dom.domRMSRes.Cl = sqrt(mean(abs(dom.domResid.Cl)).^2);
        dom.domRMSRes.Cn = sqrt(mean(abs(dom.domResid.Cn)).^2);

        dom.rmsMatrixOLS.Cx(j,i) = dom.domRMSRes.Cx;
        dom.rmsMatrixOLS.Cy(j,i) = dom.domRMSRes.Cy;
        dom.rmsMatrixOLS.Cl(j,i) = dom.domRMSRes.Cl;
        dom.rmsMatrixOLS.Cn(j,i) = dom.domRMSRes.Cn;

    end
end

for i = 1:dom.numVariables.Cz
    dom.ParamOLS.Cz = OLS.ID.Param.ThetaCz;
    dom.ParamOLS.Cm = OLS.ID.Param.ThetaCm;

    for j = 1:numDataPoints
        % OLS
        dom.ParamOLS.Cz(i) = dom.dataMatrixOLS.Cz(j, i);
        dom.ParamOLS.Cm(i) = dom.dataMatrixOLS.Cm(j, i);

        dom.domModel.Cz = OLS.ID.X.XCz_datasets * dom.ParamOLS.Cz;
        dom.domModel.Cm = OLS.ID.X.XCm_datasets * dom.ParamOLS.Cm;

        dom.domResid.Cz = OLS.ID.Y.YCz_datasets - dom.domModel.Cz;
        dom.domResid.Cm = OLS.ID.Y.YCm_datasets - dom.domModel.Cm;

        dom.domRMSRes.Cz = sqrt(mean(abs(dom.domResid.Cz)).^2);
        dom.domRMSRes.Cm = sqrt(mean(abs(dom.domResid.Cm)).^2);

        dom.rmsMatrixOLS.Cz(j,i) = dom.domRMSRes.Cz;
        dom.rmsMatrixOLS.Cm(j,i) = dom.domRMSRes.Cm;

    end
end

%%
% RMS Residual Boxplot Spread per Model Term
figure('Position', [100, 100, 700, 500]);

subplot(2,3,1)
boxplot(dom.rmsMatrixOLS.Cx, 'Positions', 1:dom.numVariables.Cx, 'Labels', cellstr(num2str((1:dom.numVariables.Cx)')));
grid on; set(gca, 'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.FontSize = 18;
set(gca,'xticklabel',{'$C_{X_0}$', '$C_{X_\alpha}$', '$C_{X_\alpha^2}$', '$C_{X_q}$', '$C_{X_{\delta_e}}$', '$C_{X_{T_c}}$'})
title("$C_X$", 'Interpreter', 'latex');
ylabel('RMS Residual');

subplot(2,3,2)
boxplot(dom.rmsMatrixOLS.Cy, 'Positions', 1:dom.numVariables.Cy, 'Labels', cellstr(num2str((1:dom.numVariables.Cy)')));
grid on; set(gca, 'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.FontSize = 18;
set(gca,'xticklabel',{'$C_{Y_0}$', '$C_{Y_\beta}$', '$C_{Y_p}$', '$C_{Y_r}$', '$C_{Y_{\delta_a}}$', '$C_{Y_{\delta_r}}$'})
title("$C_Y$", 'Interpreter', 'latex');
ylabel('RMS Residual');

subplot(2,3,3)
boxplot(dom.rmsMatrixOLS.Cz, 'Positions', 1:dom.numVariables.Cz, 'Labels', cellstr(num2str((1:dom.numVariables.Cz)')));
grid on; set(gca, 'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.FontSize = 18;
set(gca,'xticklabel',{'$C_{Z_0}$', '$C_{Z_\alpha}$', '$C_{Z_q}$', '$C_{Z_{\delta_e}}$', '$C_{Z_{T_c}}$'})
title("$C_Z$", 'Interpreter', 'latex');
ylabel('RMS Residual');

subplot(2,3,4)
boxplot(dom.rmsMatrixOLS.Cl, 'Positions', 1:dom.numVariables.Cl, 'Labels', cellstr(num2str((1:dom.numVariables.Cl)')));
grid on; set(gca, 'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.FontSize = 18;
set(gca,'xticklabel',{'$C_{l_0}$', '$C_{l_\beta}$', '$C_{l_p}$', '$C_{l_r}$', '$C_{l_{\delta_a}}$', '$C_{l_{\delta_r}}$'})
title("$C_l$", 'Interpreter', 'latex');
ylabel('RMS Residual');

subplot(2,3,5)
boxplot(dom.rmsMatrixOLS.Cm, 'Positions', 1:dom.numVariables.Cm, 'Labels', cellstr(num2str((1:dom.numVariables.Cm)')));
grid on; set(gca, 'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.FontSize = 18;
set(gca,'xticklabel',{'$C_{m_0}$', '$C_{m_\alpha}$', '$C_{m_q}$', '$C_{m_{\delta_e}}$', '$C_{m_{T_c}}$'})
title("$C_m$", 'Interpreter', 'latex');
ylabel('RMS Residual');

subplot(2,3,6)
boxplot(dom.rmsMatrixOLS.Cn, 'Positions', 1:dom.numVariables.Cn, 'Labels', cellstr(num2str((1:dom.numVariables.Cn)')));
grid on; set(gca, 'TickLabelInterpreter', 'latex');
ax = gca; ax.XAxis.FontSize = 18;
set(gca,'xticklabel',{'$C_{n_0}$', '$C_{n_\beta}$', '$C_{n_p}$', '$C_{n_r}$', '$C_{n_{\delta_a}}$', '$C_{n_{\delta_r}}$'})
title("$C_n$", 'Interpreter', 'latex');
ylabel('RMS Residual');
