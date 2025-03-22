%% ********************************************************************************************************************************
% Authors: Pietro Bianchini (001135241) Marco Valentino Caravella (001134954)
% Title of the paper: Exam Project for Structural Macroeconometrics 
% Date: 20/02/2025
% 
%% ********************************************************************************************************************************
clear all; close all; 
 
% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 12);

% Specify sheet and range
opts.Sheet = "FINAL";
opts.DataRange = "A2:L229";

% Specify column names and types
opts.VariableNames = ["DATES", "GDP", "TAX", "G", "TB3MS", "CPI_PIQ4", "MUNI1Y", "PDVMILY", "DTFP_UTIL", "HAMILTON3YP", "RESID08", "TAXNARRATIVE"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
Data = readtable("C:\Users\user\Documents\Unibo\Economics and Econometrics\Second Year\Second Semester\Structural Macroeconometrics\Matlab\Bianchini Caravella (2025)\CK_RESTUD_DATASET.xlsx", opts, "UseExcel", false);
% line 1: Dates
% line 2: GDP (in logs, per capita)
% line 3: Tax Revenue (in logs, per capita)
% line 4: Gov. Spending (in logs, per capita)
% line 5: Interest Rate
% line 6: Inflation 
% line 7: News about Tax Shocks
% line 8: News in Gov. Defense Spending
% line 9: Utilization-adjusted TFP
% line 10: Oil Price Shocks
% line 11: Monetary Policy Shocks
% line 12: Narrative Shocks

DATA = table2array(Data);

save('DATA');