%% ********************************************************************************************************************************
% Authors: Pietro Bianchini (001135241) Marco Valentino Caravella (001134954)
% Title of the paper: Exam Project for Structural Macroeconometrics 
% Date: 20/02/2025
% This code replicates Figure 1
%% ********************************************************************************************************************************
clear; close all; format bank  

load DATA

% Load your data (assuming it's already in a matrix called "data")
years = DATA(:,1);  
values = DATA(:,12); 

% Convert fractional years into datetime format
t = datetime(floor(years),1,1) + calquarters(round(4 * mod(years, 1)));


% Now t contains the actual time points corresponding to each observation.
plot(t, values, 'Color', [0 0 0.5], 'LineWidth', 1)  
xlabel('Year')
ylabel('')
title('Narrative Tax Shocks')
hold on

% Define the range of y-values
y_min = floor(min(values) / 0.1) * 0.1; 
y_max = ceil(max(values) / 0.1) * 0.1;  

% Add horizontal lines every 0.004
for y = y_min:0.1:y_max
    yline(y, 'Color', [0.8 0.8 0.8]); 
end

hold off