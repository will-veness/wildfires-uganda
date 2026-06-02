clear;
clc;
close all;
%% 

% Script to do the analysis on timeseries data from burned pixel. file
% 335cf2591a05643e298f6ac0dd6adfb2.nc
% Input data is from 1974 to 2024 to get 50 years of data

%% PIXEL 1 loading and tidy the ERA5 data 
data = readtable("pixel_lat_040_lon_2980.csv"); %time series vs total precipitaion (m)
data = rmmissing(data);
data.tp = 1000*data.tp; %as originally monthly total precip in [m], now in [mm]
pixel = sortrows(data,1,'ascend'); %ensuring in datetime order

minvals = mink(data.tp, 5);
%% SPI calculation

% Extract the precip and date part of each timestamp
SPI_monthlyprecip = table2array (pixel(:,2)); %extracting precipitation only
burned_dates = dateshift(pixel.Time, 'start', 'month');

% Group tp values by unique dates
[uniqueDates, ~, groupIndices] = unique(burned_dates);

% SPI based on precipitation following Gamma Distribution
% Remove zero precipitation values
%SPIprecipitation_nonzero = SPI_monthlyprecip(SPI_monthlyprecip > 0); Note
%this didn't change anything therefore keep SPI_monthlyprecip

% Fit gamma distribution and extract parameters
params = gamfit(SPI_monthlyprecip);
k = params(1); % Shape parameter
theta = params(2); % Scale parameter

% compute the Gamma CDF - compute cumulative probabilities for all precipitation data
F = gamcdf(SPI_monthlyprecip, k, theta);

% compute proportion of zero precip events
q = sum(SPI_monthlyprecip == 0) / length(SPI_monthlyprecip);

% Compute cumulative probabilities for nonzero precipitation
F_nonzero = gamcdf(SPI_monthlyprecip, k, theta);

% Initialize adjusted CDF array
F_adjusted = zeros(size(SPI_monthlyprecip));

% Adjust cumulative probabilities
F_adjusted(SPI_monthlyprecip > 0) = q + (1 - q) * F_nonzero;
F_adjusted(SPI_monthlyprecip == 0) = q; % Assign probability for zero precipitation

% Plot observed data and fitted Gamma CDF
x = linspace(0, max(SPI_monthlyprecip), 100); % Range of precipitation values
y = gamcdf(x, k, theta); % Fitted Gamma CDF

figure;
plot(x, y, 'r-', 'LineWidth', 2);
hold on;
scatter(SPI_monthlyprecip, F_adjusted, 'bo');
xlabel('Precipitation');
ylabel('Cumulative Probability');
legend('Gamma CDF', 'Observed Data');
grid on;

%%
spi = norminv(F_adjusted, 0, 1); % Transform to standard normal distribution

SPI_dates = table(uniqueDates, spi);

%%
figure;
plot(uniqueDates,spi);
xlabel('Date');
ylabel('SPI');
grid on;

i00 = find(SPI_dates.uniqueDates == datetime(2011,12,01));
fprintf('Monthly SPI value for December 2011 is');
disp(SPI_dates.spi(i00));
i = find(SPI_dates.uniqueDates == datetime(2012,01,01));
fprintf('Monthly SPI value for January 2012 is');
disp(SPI_dates.spi(i));
i2 = find(SPI_dates.uniqueDates == datetime(2012,02,01));
fprintf('\n Monthly SPI value for February 2012 is');
disp(SPI_dates.spi(i2));
 