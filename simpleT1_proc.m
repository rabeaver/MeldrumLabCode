clear
close all
clc

% User parameters/Load data

% Input data location and filestem, number of profiles, depths per profile,
% and echoes per depth

% Scripts folder containing T1/T2 fitting routines
addpath('/Users/tyler/Dropbox/Coding/Matlab Processing/Scripts/')

cd('/Users/tyler/Documents/Research/Aachen/Data/Museum Ludwig, Jan 2012/17.01.2012/vandongen_ana_pos1_t1/2/')
filename = 'data.dat';
year = 1909;

data = load(filename);
omitpoints = 1;

time = data(omitpoints+1:end,1);
signal = data(omitpoints+1:end,2);

% Fit to exp decay

y0_guess = 0;
A_guess = max(signal);
t1_guess = time(5); %ms

guesses = [y0_guess;A_guess;t1_guess];

CI = 90; %desired confidence interval in percent

[xfit,ypred,coeffs,residuals] = T1_fit(time,signal,guesses,CI);
% 
% figure
% subplot(2,1,1)
% hold all
% scatter(time,signal,'or')
% plot(xfit,ypred)
% subplot(2,1,2)
% hold on
% plot(time,residuals,'-b')
% plot(time,zeros(length(time)),'-k')

fid = fopen('/Users/tyler/Documents/Research/Aachen/Data/Museum Ludwig, Jan 2012/all_data.txt', 'a');
fprintf(fid, '%s\t%6f\t%6f\t%6f\t%6f\t%6f\t%6f\t%4i\t%6f\t%6f\n', filename, coeffs(1,1), coeffs(1,2), coeffs(2,1), coeffs(2,2), coeffs(3,1), coeffs(3,2),year,1/coeffs(3,1),coeffs(3,2)/(coeffs(3,1)^2));
fclose(fid);