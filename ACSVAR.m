%% ********************************************************************************************************************************
% Authors: Pietro Bianchini (001135241) Marco Valentino Caravella (001134954)
% Title of the paper: Exam Project for Structural Macroeconometrics 
% Date: 20/02/2025
% This code replicates Figure 9
%% ********************************************************************************************************************************

clc
clear
close all

global NLags TAll Sigma_Eta sig mS

%% Dataset
DataTable = readtable('CK_RESTUD_DATASET.xlsx', 'Sheet', 4);
DataSet = [DataTable{:,2}, DataTable{:,3}, DataTable{:,4}, DataTable{:,12}];
AllDataSet = DataSet;
M = size(DataSet,2); 

%% VAR
NLags = 4;
TAll = size(DataSet,1) - NLags;
BootstrapIterations = 1000;
options = optimset('MaxFunEvals',200000, 'TolFun',1e-500, 'MaxIter',200000, 'TolX',1e-50);

%% Model dimensions for Partial Shock Identification
% n: variabili endogene (3: GDP, TAX, G)
% r: dimensione del proxy (1: TAXNARRATIVE)
r = 1;
n = M - r;  % 3
g_shock = 1;

%% Useful Matrices
DuplicationMatrix = DuplicationMatrixFunction(M);
mDD = (DuplicationMatrix' * DuplicationMatrix)^(-1) * DuplicationMatrix';
mNN = DuplicationMatrix * mDD;

KommutationMatrix = CommutationMatrixFunction(M);  % Dimensioni: 16x16 se M=4
NMatrix = 0.5 * (eye(M^2) + KommutationMatrix);

DuplicationMatrixr = DuplicationMatrixFunction(r);
mDDg = (DuplicationMatrixr' * DuplicationMatrixr)^(-1) * DuplicationMatrixr;

DuplicationMatrixN = DuplicationMatrixFunction(n);
mDDn = (DuplicationMatrixN' * DuplicationMatrixN)^(-1) * DuplicationMatrixN';

% KommutationMatrixNg: dimensione (n*g_shock) x (n*g_shock) = 3x3
KommutationMatrixNg = zeros(n*g_shock, n*g_shock);
KommutationMatrixNg(1,1) = 1;
KommutationMatrixNg(2,3) = 1;
KommutationMatrixNg(3,2) = 1;

%% Exploratory data analysis

numObs = size(DataSet,1);
dates = datetime(1950,1,1) + calquarters(0:numObs-1);

% Plot time series
figure('Name','Figure 1: Main Variables','Units','normalized','Position',[0.1 0.1 0.8 0.8]);
subplot(3,1,1)
plot(dates, AllDataSet(:,1), 'LineWidth',1.2);
title('GDP');
xlabel('Year');
ylabel('GDP');
xtickformat('yyyy');
grid on;

subplot(3,1,2)
plot(dates, AllDataSet(:,2), 'LineWidth',1.2);
title('TAX');
xlabel('Year');
ylabel('Tax');
xtickformat('yyyy');
grid on;

subplot(3,1,3)
plot(dates, AllDataSet(:,3), 'LineWidth',1.2);
title('Government Spending (G)');
xlabel('Year');
ylabel('G');
xtickformat('yyyy');
grid on;

figure('Name','Figure 2: Tax Narrative','Units','normalized','Position',[0.3 0.3 0.4 0.4]);
plot(dates, AllDataSet(:,4), 'LineWidth',1.2);
title('Tax Narrative');
xlabel('Year');
ylabel('Tax Narrative');
xtickformat('yyyy');
grid on;

%ADF test su (GDP, TAX, G)
alpha = 0.05;
nVars = 3;
varNamesEDA = {'GDP','TAX','G'};
hADF = zeros(nVars,1);
pADF = zeros(nVars,1);
decisionADF = cell(nVars,1);

for i = 1:nVars
    [h, pValue] = adftest(AllDataSet(:,i));
    hADF(i) = h;
    pADF(i) = pValue;
    if pValue < alpha
        decisionADF{i} = 'Rejected';  % stazionaria
    else
        decisionADF{i} = 'Not Rejected';  % Non stazionaria
    end
end

T_ADF = table(varNamesEDA', pADF, decisionADF, 'VariableNames', {'Variable','pValue','ADF_Test'});
disp('--- ADF Test Results ---');
disp(T_ADF);
writetable(T_ADF, 'ADF_Results_OurDataset.csv');

%% VAR Estimation
VAR_Const = [NaN; NaN; NaN];
VAR_A = {NaN(3,3), NaN(3,3), NaN(3,3), NaN(3,3)};

Mdl = varm('Constant', VAR_Const, 'AR', VAR_A);
[EstMdl, ~, logL_Rest, Errors] = estimate(Mdl, AllDataSet(:,1:3));

A_Const = EstMdl.Constant;  % 3x1
A_1 = EstMdl.AR{1};         % 3x3
A_2 = EstMdl.AR{2};         % 3x3
A_3 = EstMdl.AR{3};         % 3x3
A_4 = EstMdl.AR{4};         % 3x3
Sigma_Eta = EstMdl.Covariance;  % 3x3
Sigma_Eta_Sample = Sigma_Eta;

%% Calcolo di Sigma_Eta_Sample_full (4x4) dall'intero dataset
Sigma_Eta_Sample_full = cov(AllDataSet);
Omega_eta_full = 2 * mDD * kron(Sigma_Eta_Sample_full, Sigma_Eta_Sample_full) * mDD';
temp = mDD * kron(Sigma_Eta_Sample_full, Sigma_Eta_Sample_full) * mDD';
temp_scaled = (1/TAll) * temp;
StandardErrors_Omega = sqrt(diag(temp_scaled));

%% Partial Shocks Identification Strategy
Sigma_u = Sigma_Eta_Sample_full(1:3, 1:3);
Sigma_uv = Sigma_Eta_Sample_full(1:3, end);
Sigma_vu = Sigma_Eta_Sample_full(end, 1:3);
Sigma_v = Sigma_Eta_Sample_full(end, end);

DiffZero = tril(Sigma_Eta_Sample_full);
VechSigmaU = DiffZero(DiffZero ~= 0);
P_matrix = eye(10);
lambda = P_matrix * VechSigmaU;
Omega_lambda = P_matrix * Omega_eta_full * P_matrix';

% Definisce il vettore dei momenti "sig" (4x1):
%   La prima componente è (Sigma_vu/Sigma_u)*Sigma_uv (uno scalare),
%   Le successive 3 componenti sono la trasposta di Sigma_vu (3x1)
sig = [ (Sigma_vu / Sigma_u) * Sigma_uv;
         Sigma_vu' ];

% Calcolo di F_lambda per il CMD (dimensioni per r = 1, n = 3)
F_lambda = [ -mDDg * kron(Sigma_vu / Sigma_u, Sigma_vu / Sigma_u) * (mDDn)', 2*mDDg * kron(Sigma_vu / Sigma_u, eye(r)), zeros(0.5*r*(r+1), 0.5*r*(r+1));
             zeros(n*r, 0.5*n*(n+1)), eye(n*r), zeros(n*r, 0.5*r*(r+1)) ];
Omega_zeta = F_lambda * Omega_lambda * F_lambda';
mS = Omega_zeta;  % Matrice dei pesi per il CMD (4x4)

%% Some tests
% AIC e BIC 
maxLag = 12;
AIC_vals = zeros(maxLag,1);
BIC_vals = zeros(maxLag,1);
for lag = 1:maxLag
    mdl_temp = varm('Constant', repmat(NaN, n, 1), 'AR', repmat({NaN(n,n)},1,lag));
    [estTemp,~,logL_temp] = estimate(mdl_temp, AllDataSet(:,1:3));
    numParams = lag * n^2 + n;  % stima approssimativa del numero di parametri
    AIC_vals(lag) = -2*logL_temp + 2*numParams;
    BIC_vals(lag) = -2*logL_temp + numParams * log(TAll);
end

% Table AIC BIC
T_IC = table((1:maxLag)', AIC_vals, BIC_vals, 'VariableNames', {'LagOrder','AIC','BIC'});
writetable(T_IC, 'InfoCriteria.csv');  % Esporta la tabella in CSV

% Jarque-Bera, Ljung-Box and ARCH
JB_p = zeros(n,1);
LB_p = zeros(n,1);
ARCH_p = zeros(n,1);
for i = 1:n
    [~, JB_p(i)] = jbtest(Errors(:,i));
    [~, LB_p(i)] = lbqtest(Errors(:,i), 'Lags', 4);
    [~, ARCH_p(i)] = archtest(Errors(:,i), 'Lags', 4);
end
T_diag = table({'GDP'; 'TAX'; 'GovSpending'}, JB_p, LB_p, ARCH_p, 'VariableNames', {'Equation','JB_p','LB_p','ARCH_p'});
writetable(T_diag, 'DiagnosticTests.csv'); 


%% CMD Estimation
% vg = [beta1 (3x1); phi] (4x1)
vg_initial = [0.7; 0; 0; 0];
[Param, Distance, exitflag, output, grad, Hessian_Matrix] = fminunc(@MinimumDistance, vg_initial, options);

phi_est = Param(4);
B_est = Param(1:3);

% Calcola SE_est da VAR_est:
F_theta_est = [ zeros(0.5*r*(r+1), n*g_shock), 2*mDDg * kron(phi_est, eye(r));
                kron(eye(n), phi_est) * KommutationMatrixNg, kron(B_est, eye(r)) ];
VAR_est = inv(F_theta_est' / mS * F_theta_est) / TAll;
SE_est = sqrt(diag(VAR_est));

%% CMD Restricted
vg_initial_restricted = [Param(1); Param(4)];  % 2x1
[Restricted_Param, Restricted_Distance, exitflag, output, grad, Restricted_Hessian_Matrix] = fminunc(@MinimumDistance_Restricted, vg_initial_restricted, options);
Partial_ExogTest = 1 - chi2cdf(Restricted_Distance*TAll - Distance*TAll, 2);

%% Bootstrap
DataSet_Sample = DataSet;
wbGS = waitbar(0, 'Running the bootstrap');
for boot = 1:BootstrapIterations
    waitbar(boot/BootstrapIterations);
    TBoot = datasample(1:TAll, TAll);
    Errors_Boot = Errors(TBoot, :);
    StoreErrors_Boot(:, :, boot) = Errors_Boot;
    
    T = TAll;
    DataSet_Bootstrap = zeros(T+NLags, n);
    DataSet_Bootstrap(1:NLags, :) = DataSet_Sample(1:NLags, 1:3);
    for t = NLags+1:T+NLags
        Y = A_Const + A_1 * DataSet_Bootstrap(t-1, :)' + ...
            A_2 * DataSet_Bootstrap(t-2, :)' + ...
            A_3 * DataSet_Bootstrap(t-3, :)' + ...
            A_4 * DataSet_Bootstrap(t-4, :)' + Errors_Boot(t-NLags, 1:3)';
        DataSet_Bootstrap(t, :) = Y';
    end
    
    DataSet_VAR = DataSet_Bootstrap;
    EstMdl_Boot = estimate(Mdl, DataSet_VAR);
    A_Const_Boot(:, :, boot) = EstMdl_Boot.Constant;
    A_1_Boot(:, :, boot) = EstMdl_Boot.AR{1};
    A_2_Boot(:, :, boot) = EstMdl_Boot.AR{2};
    A_3_Boot(:, :, boot) = EstMdl_Boot.AR{3};
    A_4_Boot(:, :, boot) = EstMdl_Boot.AR{4};
    Sigma_Eta = EstMdl_Boot.Covariance;
    Sigma_Eta_Boot(:, :, boot) = Sigma_Eta;
    
    InitialValue = Param;
    [Param_Boot, ~, ~, ~, ~, ~] = fminunc(@MinimumDistance, InitialValue, options);
    Param_Boot_All(:, boot) = Param_Boot;
    
    G_REST_Boot(:, :, boot) = Param_Boot(1:3);  % Impatto dello shock come beta1 (3x1)
end
close(wbGS);

%% IRF Calculation
HorizonIRF = 24;
J = [eye(n) zeros(n, n*(NLags-1))];
CompanionMatrix = [A_1 A_2 A_3 A_4;
                   eye(n*(NLags-1)) zeros(n*(NLags-1), n)];
IRF_Sample = zeros(n, 1, HorizonIRF+1);
for h = 0:HorizonIRF
    IRF_Sample(:, 1, h+1) = J * (CompanionMatrix^h) * J' * G_REST_Boot(:, :, 1);
end

IRF_Boot = zeros(n, 1, HorizonIRF+1, BootstrapIterations);
for i = 1:BootstrapIterations
    CompMatrix_Boot = [A_1_Boot(:, :, i) A_2_Boot(:, :, i) A_3_Boot(:, :, i) A_4_Boot(:, :, i)];
    CompanionMatrix_Boot = [CompMatrix_Boot;
                           eye(n*(NLags-1)) zeros(n*(NLags-1), n)];
    for h = 0:HorizonIRF
        IRF_Boot(:, 1, h+1, i) = J * (CompanionMatrix_Boot^h) * J' * G_REST_Boot(:, :, i);
    end
end

%% Plot IRFs
LineWidth_IRF = 2;
FontSizeIRFGraph = 12;
alpha_coverage = 10;
figure(1)
Titles = {'GDP Response','Tax Response','Government Spending Response'};
for i = 1:n
    IRF_i = squeeze(IRF_Sample(i, 1, :));
    IRF_Boot_i = squeeze(IRF_Boot(i, 1, :, :));  % (HorizonIRF+1) x BootstrapIterations
    SE_i = sqrt(mean(IRF_Boot_i.^2, 2) - mean(IRF_Boot_i, 2).^2);
    for j = 1:BootstrapIterations
        m_max(j,:) = max(abs(IRF_Boot_i(:, j) - IRF_i) ./ SE_i);
    end
    lower_bound = IRF_i - SE_i .* prctile(m_max, 100 - alpha_coverage);
    upper_bound = IRF_i + SE_i .* prctile(m_max, 100 - alpha_coverage);
    
    subplot(n,1,i)
    x = 0:HorizonIRF;
    fill([x, fliplr(x)], [lower_bound', fliplr(upper_bound')], [0 0.4470 0.7410], 'linestyle', 'none', 'facealpha', 0.3);
    hold on
    plot(x, IRF_i, 'Color', [0 0.4470 0.7410], 'LineWidth', LineWidth_IRF);
    plot(x, zeros(size(x)), 'k', 'LineWidth', 1);
    title(Titles{i}, 'FontSize', FontSizeIRFGraph);
    ylabel('Response','FontSize',FontSizeIRFGraph);
    xlabel('Horizon (quarters)','FontSize',FontSizeIRFGraph);
    axis tight
end

%% Stima di G1
paramNames = {'Beta1_1'; 'Beta1_2'; 'Beta1_3'; 'Phi'};
Estimates = num2cell(Param);
StdErrors = num2cell(SE_est);  % SE_est deve essere già definito
T_G1 = table(paramNames, Estimates, StdErrors, 'VariableNames', {'Parameter', 'Estimate', 'StdError'});
writetable(T_G1, 'EstimatedG1.csv');  % Esporta in CSV per poi usarlo in LaTeX

%% Plot individual IRFs
for i = 1:n
    figure;
    IRF_i = squeeze(IRF_Sample(i, 1, :));  % Vettore IRF per la variabile i
    IRF_Boot_i = squeeze(IRF_Boot(i, 1, :, :));  % (HorizonIRF+1) x BootstrapIterations
    SE_i = sqrt(mean(IRF_Boot_i.^2, 2) - mean(IRF_Boot_i, 2).^2);
    for j = 1:BootstrapIterations
        m_max(j,:) = max(abs(IRF_Boot_i(:, j) - IRF_i) ./ SE_i);
    end
    lower_bound = IRF_i - SE_i .* prctile(m_max, 100 - alpha_coverage);
    upper_bound = IRF_i + SE_i .* prctile(m_max, 100 - alpha_coverage);
    
    x = 0:HorizonIRF;
    fill([x, fliplr(x)], [lower_bound', fliplr(upper_bound')], [0 0.4470 0.7410], 'linestyle', 'none', 'facealpha', 0.3);
    hold on
    plot(x, IRF_i, 'Color', [0 0.4470 0.7410], 'LineWidth', LineWidth_IRF);
    plot(x, zeros(size(x)), 'k', 'LineWidth', 1);
    title(sprintf('IRF of %s to a Tax Shock', getVariableName(i)), 'FontSize', FontSizeIRFGraph);
    ylabel('Response', 'FontSize', FontSizeIRFGraph);
    xlabel('Horizon (quarters)', 'FontSize', FontSizeIRFGraph);
    axis tight;
    
    % Salva la figura in formato PDF per LaTeX (opzionale)
    print(gcf, sprintf('IRF_Var%d.pdf', i), '-dpdf', '-r300');
end

function varName = getVariableName(idx)
    switch idx
        case 1
            varName = 'GDP';
        case 2
            varName = 'TAX';
        case 3
            varName = 'Government Spending';
        otherwise
            varName = 'Unknown';
    end
end

%% Display Results
clc
disp('--- ADF Test Results ---');
disp(T_ADF);
disp('************************************************************')
disp('Partial Shocks Identification Results')
disp('----------------------------------------')
disp('Estimated parameters (Param):')
disp(Param')
disp('Standard Errors:')
disp(SE_est)
disp('Partial Exogeneity Test Statistic:')
test_stat = Restricted_Distance * TAll - Distance * TAll;
disp(test_stat)
disp('p-value:')
disp(1 - chi2cdf(test_stat, 2))
disp('----------------------------------------')

% Informative Criteria (AIC e BIC)
if isfile('InfoCriteria.csv')
    T_IC = readtable('InfoCriteria.csv');
    disp('Information Criteria (AIC, BIC) for different lag orders:')
    disp(T_IC)
else
    disp('InfoCriteria.csv not found.')
end

% Diagnostic Tests sui residui del VAR
if isfile('DiagnosticTests.csv')
    T_diag = readtable('DiagnosticTests.csv');
    disp('Diagnostic Tests on VAR Residuals (p-values):')
    disp(T_diag)
else
    disp('DiagnosticTests.csv not found.')
end

% Risultati della stima di G1 e relativi errori standard
if isfile('EstimatedG1.csv')
    T_G1 = readtable('EstimatedG1.csv');
    disp('Estimated G1 parameters and Standard Errors:')
    disp(T_G1)
else
    disp('EstimatedG1.csv not found.')
end

disp('************************************************************')


%% Local Functions 

function K = CommutationMatrixFunction(n)
    % Costruisce la matrice di commutazione K di dimensione n^2 x n^2
    K = zeros(n*n, n*n);
    for i = 1:n
        for j = 1:n
            index1 = (j-1)*n + i;
            index2 = (i-1)*n + j;
            K(index2, index1) = 1;
        end
    end
end

function D = DuplicationMatrixFunction(n)
    % Genera la matrice di duplicazione per matrici simmetriche n x n.
    vechindex = @(r,c,n) n*r + c - r*(r+1)/2;
    D = zeros(n*n, n*(n+1)/2);
    for j = 1:n*n
        rowo = floor((j-1)/n);
        colo = mod(j-1, n);
        if (rowo <= colo)
            D(j, vechindex(rowo, colo, n)+1) = 1;
        else
            D(j, vechindex(colo, rowo, n)+1) = 1;
        end
    end
end

function yv = vec(Y)
    [n, m] = size(Y);
    yv = reshape(Y, n*m, 1);
end

function v = vech(matA)
    [M, N] = size(matA);
    if M ~= N
        error('Input must be a symmetric matrix.');
    end
    v = [];
    for ii = 1:M
        v = [v; matA(ii:end, ii)];
    end
end

function [MinDist] = MinimumDistance(vg)
    % vg è un vettore 4x1: [beta1 (3x1); phi]
    global mS
    global sig  % sig deve essere un vettore 4x1
    beta1 = vg(1:3);  % 3x1
    phi = vg(4);      % scalare
    ftheta = [phi^2; phi * beta1];  % 4x1
    diff = sig - ftheta;
    MinDist = diff' * inv(mS) * diff;
end

function [MinDist] = MinimumDistance_Restricted(vg)
    % vg è un vettore 2x1: [beta1_restricted; phi]
    global mS
    global sig  % per la restrizione consideriamo le prime 2 componenti di sig
    beta1 = vg(1);  % scalare
    phi = vg(2);    % scalare
    ftheta = [phi^2; phi * beta1];  % 2x1
    diff = sig(1:2) - ftheta;
    MinDist = diff' * inv(mS(1:2,1:2)) * diff;
end
