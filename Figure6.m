%% ********************************************************************************************************************************
% Authors: Pietro Bianchini (001135241) Marco Valentino Caravella (001134954)
% Title of the paper: Exam Project for Structural Macroeconometrics 
% Date: 20/02/2025
% This code replicates Figure 6
%% ********************************************************************************************************************************
clear; close all;addpath('additional files'); format bank  

load DATA

% Repositioning of the variables
DATA = DATA(:, [1, 3, 4, 2, 5, 6, 7, 8, 9, 10, 11, 12]);
Data = movevars(Data, 'GDP', 'After', Data.Properties.VariableNames{4});

%  Proxy SVAR with Fed Funds Rates and Inflation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 dates = DATA(:,1);
 VARNA.vars  = [DATA(:,2:6)];
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
  VARNA = doPVAR2(VARNA);
  nboot = 10000; % Number of replications
  clevel = 95; % Confidence Level
  VARNAbs = doPVARbs2(VARNA,DAT,nboot,clevel);

  % Parameter Estimates
  fprintf('thetaG = %f \n',VARNA.thetaG);
  fprintf('thetaY = %f \n',VARNA.thetaY);
  fprintf('gammaT = %f \n',VARNA.gammaT);
  fprintf('zetaT = %f \n',VARNA.zetaT);
  fprintf('zetaG = %f \n',VARNA.zetaG);

  % Plot Impulse Response to a Tax Shock
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FIG.axes = [-2 5.5;-10 10;-10 10];

  f=figure;
    
    box on
        plot(VARNA.irs(:,3),'LineWidth',2,'Color', [0 0 0.5] );
        hold on
        plot(VARNAbs.irsH(:,3),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        plot(VARNAbs.irsL(:,3),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        axis([1 20 FIG.axes(1,1) FIG.axes(1,2)])
        hline(0,'k-')
        ti=title('GDP');
        xl=xlabel('quarters');
        yl=ylabel('percent');
    
        set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
        set([ti], 'FontName', 'AvantGarde','FontSize',16); 
 
  saveas(gcf, '6GDP.png');         
   
    f=figure;
    
    box on
        plot(VARNA.irs(:,2),'LineWidth',2,'Color', [0 0 0.5]);
        hold on       
        plot(VARNAbs.irsH(:,2),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        plot(VARNAbs.irsL(:,2),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        axis([1 20 FIG.axes(2,1) FIG.axes(2,2)])
        hline(0,'k-')
        ti=title('Government Spending');
        xl=xlabel('quarters');
        yl=ylabel('percent');
    
        set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
        set([ti], 'FontName', 'AvantGarde','FontSize',16); 
         
   saveas(gcf, '6GovernmentSpending.png'); 

   f=figure;
    
    box on
        plot(VARNA.irs(:,1),'LineWidth',2,'Color', [0 0 0.5]);
        hold on       
        plot(VARNAbs.irsH(:,1),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        plot(VARNAbs.irsL(:,1),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        axis([1 20 FIG.axes(3,1) FIG.axes(3,2)])
        hline(0,'k-')
        ti=title('Tax Revenues');
        xl=xlabel('quarters');
        yl=ylabel('percent');
    
        set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
        set([ti], 'FontName', 'AvantGarde','FontSize',16); 
     
    saveas(gcf, '6TaxRevenues.png');

    f=figure;
    
    box on
        plot(VARNA.irs(:,4),'LineWidth',2,'Color', [0 0 0.5]);
        hold on       
        plot(VARNAbs.irsH(:,4),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        plot(VARNAbs.irsL(:,4),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        axis([1 20 FIG.axes(3,1) FIG.axes(3,2)])
        hline(0,'k-')
        ti=title('Interest Rate');
        xl=xlabel('quarters');
        yl=ylabel('percent');
    
        set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
        set([ti], 'FontName', 'AvantGarde','FontSize',16); 

saveas(gcf, '6InterestRate.png'); 

   f=figure;
    
    box on
        plot(VARNA.irs(:,5),'LineWidth',2,'Color', [0 0 0.5]);
        hold on       
        plot(VARNAbs.irsH(:,5),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        plot(VARNAbs.irsL(:,5),'LineWidth',1,'Color',[0 0 0.5],'LineStyle','--' ); 
        hold on
        axis([1 20 FIG.axes(3,1) FIG.axes(3,2)])
        hline(0,'k-')
        ti=title('Inflation');
        xl=xlabel('quarters');
        yl=ylabel('percent');
    
        set([xl,yl], 'FontName', 'AvantGarde','FontSize',14);
        set([ti], 'FontName', 'AvantGarde','FontSize',16); 

   saveas(gcf, '6Inflation.png'); 