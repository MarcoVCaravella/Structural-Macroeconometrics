%% ********************************************************************************************************************************
% Authors: Pietro Bianchini (001135241) Marco Valentino Caravella (001134954)
% Title of the paper: Exam Project for Structural Macroeconometrics 
% Date: 20/02/2025
% This code replicates Figure 5
%% ********************************************************************************************************************************
clear; close all;addpath('additional files'); format bank  

load DATA

% Repositioning of the variables
DATA = DATA(13:end, [1, 3, 4, 2, 5, 6, 7, 8, 9, 10, 11, 12]);
Data = movevars(Data, 'GDP', 'After', Data.Properties.VariableNames{4});

%  Proxy SVAR - Benchmark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dates = DATA(:,1);
VARNA.vars  = [DATA(:,2:4)];
VARNA.p = 4;
[T,n]  = size(VARNA.vars);
VARNA.DET = [ones(T,1) (1:T)' ((1:T).^2)' dates==1975.25];
VARNA.irhor = 20;

VARNA.taxshocks = [DATA(:,12)];
 
% Demean the narrative shocks
 for j=1:size(VARNA.taxshocks,2)
    VARNA.taxshocks(VARNA.taxshocks(:,j)~=0,j)=(VARNA.taxshocks(VARNA.taxshocks(:,j)~=0,j)-mean(VARNA.taxshocks(VARNA.taxshocks(:,j)~=0,j)));
 end
  
DAT.TRY = 0.182160825; % Average ratio of federal tax revenues to GDP
DAT.GY  = 0.204839843; % Average ratio of federal expenditures to GDP

VARNA.tshocksize =  0.01/DAT.TRY;
VARNA.gshocksize =  0.01/DAT.GY; 
VARNA = doPVAR(VARNA);
nboot = 10000; % Number of replications
clevel = 95; % Confidence Level
VARNAbs = doPVARbs(VARNA,DAT,nboot,clevel);

% Parameter Estimates
fprintf('thetaG = %f \n',VARNA.thetaG);
fprintf('thetaY = %f \n',VARNA.thetaY);
fprintf('gammaT = %f \n',VARNA.gammaT);
fprintf('zetaT = %f \n',VARNA.zetaT);
fprintf('zetaG = %f \n',VARNA.zetaG);

% Store results from the first case
irs_case1 = VARNA.irs;
irsH_case1 = VARNAbs.irsH;
irsL_case1 = VARNAbs.irsL;

%  Proxy SVAR with MUNI1Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dates = DATA(:,1);
VARNA.vars  = [DATA(:,2:4)];
VARNA.p = 4;
[T,n]  = size(VARNA.vars);
VARNA.DET = [ones(T,1) (1:T)' ((1:T).^2)' dates==1975.25 DATA(:,7)];
VARNA.irhor = 20;

VARNA.taxshocks = [DATA(:,12)];
 
% Demean the narrative shocks
 for j=1:size(VARNA.taxshocks,2)
    VARNA.taxshocks(VARNA.taxshocks(:,j)~=0,j)=(VARNA.taxshocks(VARNA.taxshocks(:,j)~=0,j)-mean(VARNA.taxshocks(VARNA.taxshocks(:,j)~=0,j)));
 end
  
DAT.TRY = 0.182160825; % Average ratio of federal tax revenues to GDP
DAT.GY  = 0.204839843; % Average ratio of federal expenditures to GDP

VARNA.tshocksize =  0.01/DAT.TRY;
VARNA.gshocksize =  0.01/DAT.GY; 
VARNA = doPVAR(VARNA);
nboot = 10000; % Note:  the paper uses 10000
clevel = 95; % Confidence Level
VARNAbs = doPVARbs(VARNA,DAT,nboot,clevel);

% Parameter Estimates
fprintf('thetaG = %f \n',VARNA.thetaG);
fprintf('thetaY = %f \n',VARNA.thetaY);
fprintf('gammaT = %f \n',VARNA.gammaT);
fprintf('zetaT = %f \n',VARNA.zetaT);
fprintf('zetaG = %f \n',VARNA.zetaG);

% Store results from the first case
irs_case2 = VARNA.irs;
irsH_case2 = VARNAbs.irsH;
irsL_case2 = VARNAbs.irsL;

%  Proxy SVAR with PDVMILY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dates = DATA(:,1);
VARNA.vars  = [DATA(:,2:4)];
VARNA.p = 4;
[T,n]  = size(VARNA.vars);
VARNA.DET = [ones(T,1) (1:T)' ((1:T).^2)' dates==1975.25 DATA(:,8)];
VARNA.irhor = 20;

VARNA.taxshocks = [DATA(:,12)];
 
% Demean the narrative shocks
 for j=1:size(VARNA.taxshocks,2)
    VARNA.taxshocks(VARNA.taxshocks(:,j)~=0,j)=(VARNA.taxshocks(VARNA.taxshocks(:,j)~=0,j)-mean(VARNA.taxshocks(VARNA.taxshocks(:,j)~=0,j)));
 end
  
DAT.TRY = 0.182160825; % Average ratio of federal tax revenues to GDP
DAT.GY  = 0.204839843; % Average ratio of federal expenditures to GDP

VARNA.tshocksize =  0.01/DAT.TRY;
VARNA.gshocksize =  0.01/DAT.GY; 
VARNA = doPVAR(VARNA);
nboot = 10000; % Note:  the paper uses 10000
clevel = 95; % Confidence Level
VARNAbs = doPVARbs(VARNA,DAT,nboot,clevel);

% Parameter Estimates
fprintf('thetaG = %f \n',VARNA.thetaG);
fprintf('thetaY = %f \n',VARNA.thetaY);
fprintf('gammaT = %f \n',VARNA.gammaT);
fprintf('zetaT = %f \n',VARNA.zetaT);
fprintf('zetaG = %f \n',VARNA.zetaG);

% Store results from the second case
irs_case3 = VARNA.irs;
irsH_case3 = VARNAbs.irsH;
irsL_case3 = VARNAbs.irsL;

% Plot Impulse Response to a Tax Shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIG.axes = [-2 5.5; -10 10; -10 10];

figure;
box on
h1 = plot(irs_case1(:,3), 'LineWidth', 2, 'Color', [0 0 0.5], 'DisplayName', 'Case 1'); hold on
plot(irsH_case1(:,3), 'LineWidth', 1, 'Color', [0 0 0.5], 'LineStyle', '--'); 
plot(irsL_case1(:,3), 'LineWidth', 1, 'Color', [0 0 0.5], 'LineStyle', '--'); 

h2 = plot(irs_case2(:,3), 'LineWidth', 2, 'Color', [0.8 0 0], 'DisplayName', 'Case 2'); hold on
plot(irsH_case2(:,3), 'LineWidth', 1, 'Color', [0.8 0 0], 'LineStyle', '--'); 
plot(irsL_case2(:,3), 'LineWidth', 1, 'Color', [0.8 0 0], 'LineStyle', '--'); 

h3 = plot(irs_case3(:,3), 'LineWidth', 2, 'Color', [0 1 0], 'DisplayName', 'Case 1'); hold on
plot(irsH_case3(:,3), 'LineWidth', 1, 'Color', [0 1 0], 'LineStyle', '--'); 
plot(irsL_case3(:,3), 'LineWidth', 1, 'Color', [0 1 0], 'LineStyle', '--'); 

axis([1 20 FIG.axes(1,1) FIG.axes(1,2)]);
hline(0,'k-');
title('GDP');
xlabel('Quarters');
ylabel('Percent');
legend([h1, h2, h3], {'Benchmark', 'MUNI1Y', 'PDVMILY'}); 
set(gca, 'FontName', 'AvantGarde', 'FontSize', 14);

saveas(gcf, '5GDP.png'); 