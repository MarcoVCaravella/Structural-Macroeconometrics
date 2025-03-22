%% ********************************************************************************************************************************
% Authors: Pietro Bianchini (001135241) Marco Valentino Caravella (001134954)
% Title of the paper: Exam Project for Structural Macroeconometrics 
% Date: 20/02/2025
% This code replicates Figure 2
%% ********************************************************************************************************************************
clear all; close all; addpath('additional files'); format bank  

load DATA
% line 1: Dates
% line 2: GDP (in logs, per capita)
% line 3: Tax Revenues (in logs, per capita)
% line 4: Government Spending (in logs, per capita)

% Repositioning of the variables
DATA = DATA(:, [1, 3, 4, 2, 5, 6, 7, 8, 9, 10, 11, 12]);
Data = movevars(Data, 'GDP', 'After', Data.Properties.VariableNames{4});

%  Blanchard Perotti SVAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 dates = DATA(:,1);
 VAR.vars = [DATA(:,2:4)];
 VAR.p =  4;
 [T,n]  = size(VAR.vars);
 VAR.DET = [ones(T,1) (1:T)' ((1:T).^2)' dates==1975.25];
 VAR.irhor = 20;
 VAR.BP.thetaY=2.08; 
 VAR.BP.gammaT=0; VAR.BP.gammaY=0;
 VAR.d0 = [-0.0382;-0.0196;0.0691;0.0217;0.0232;0.0086];
 DAT.TRY = 0.182160825; % Average ratio of federal tax revenues to GDP
 DAT.GY  = 0.204839843; % Average ratio of federal expenditures to GDP
 VAR.tshocksize =  0.01/DAT.TRY;
 VAR.gshocksize =  0.01/DAT.GY; 
 
 VAR = doVAR(VAR);
 nboot = 10000; % Number of replications
 clevel = 95; % Confidence Level
 VARbs = doVARbs(VAR,nboot,clevel);
  
 % Parameter Estimates
 fprintf('thetaG = %f \n',VAR.thetaG);
 fprintf('thetaY = %f \n',VAR.thetaY);
 fprintf('gammaT = %f \n',0);
 fprintf('zetaT = %f \n',VAR.zetaT);
 fprintf('zetaG = %f \n',VAR.zetaG);

 % Plot Impulse Response to a Tax Shock
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 FIG.axes = [-2 5.5;-10 10;-10 10];

 f=figure;
    
    box on
        plot(VAR.irs(:,3),'LineWidth',2,'Color', [0 0 0.5] );
        hold on
        plot(VARbs.irsH(:,3),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        plot(VARbs.irsL(:,3),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        axis([1 20 FIG.axes(1,1) FIG.axes(1,2)])
        hline(0,'k-')
        ti=title('Output');
        xl=xlabel('quarters');
        yl=ylabel('percent');
    
        set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
        set([ti], 'FontName', 'AvantGarde','FontSize',16); 

saveas(gcf, '2GDP.png'); 

f=figure;
    
    box on
        plot(VAR.irsg(:,3),'LineWidth',2,'Color', [0 0 0.5] );
        hold on
        plot(VARbs.irsgH(:,3),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        plot(VARbs.irsgL(:,3),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        axis([1 20 FIG.axes(1,1) FIG.axes(1,2)])
        hline(0,'k-')
        ti=title('Output');
        xl=xlabel('quarters');
        yl=ylabel('percent');
    
        set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
        set([ti], 'FontName', 'AvantGarde','FontSize',16);