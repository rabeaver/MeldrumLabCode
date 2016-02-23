clear
close all
clc

%% User parameters/Load data

% Input data location and filestem, number of profiles, depths per profile,
% and echoes per depth

% Scripts folder containing T1/T2 fitting routines
% addpath('/Users/tyler/Dropbox/Coding/Matlab Processing/Scripts/')

datadir = '/Users/tyler/Desktop/';
parfilestem = 'acqu';
datafilestem = 'data';

params.acqTime = readpar_Kea(strcat(datadir,parfilestem,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(datadir,parfilestem,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(datadir,parfilestem,'.par'),'nrScans');
params.rxPhase = readpar_Kea(strcat(datadir,parfilestem,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(datadir,parfilestem,'.par'),'rxGain');
params.nrPtsT1 = readpar_Kea(strcat(datadir,parfilestem,'.par'),'nrPntsT1');
params.t1Max = readpar_Kea(strcat(datadir,parfilestem,'.par'),'tMax');
params.t1Est = readpar_Kea(strcat(datadir,parfilestem,'.par'),'t1Est');
params.repTime = readpar_Kea(strcat(datadir,parfilestem,'.par'),'repTime');
params.b1Freq = readpar_Kea(strcat(datadir,parfilestem,'.par'),'b1Freq');

% Make T1 time vector
amax = 1;
amin = exp(-params.t1Max/params.t1Est);
astep = (amax-amin)/(params.nrPtsT1-1);
t1TimeAxis = zeros(params.nrPtsT1,1);

for i = 1:1:params.nrPtsT1
    t1TimeAxis(i) = params.t1Est*(-log(amax-astep*(i-1)));
end

% Load data
amplitude = dlmread(strcat(datadir,datafilestem,'.dat'));
amplitude = amplitude(:,2)./max(amplitude(:,2));

%% Fit to exp decay

y0_guess = 0;
A_guess = 1;
t1_guess = 20; %ms

guesses = [y0_guess;A_guess;t1_guess];

CI = 90; %desired confidence interval in percent

[xfit,ypred,coeffs,residuals,J,mse,ci,se] = T1fit(t1TimeAxis(1:end),amplitude(1:end),guesses,CI);

figure
subplot(2,1,1)
hold all
scatter(t1TimeAxis,amplitude,'or')
plot(xfit,ypred)
subplot(2,1,2)
hold on
plot(t1TimeAxis,residuals,'-b')
plot(t1TimeAxis,zeros(length(t1TimeAxis)),'-k')


% sprintf('A = %f +/- %f',coeffs(1),coeffs_err(1))
% sprintf('T2 = %f +/- %f',coeffs(3),ci(3))
sprintf('T2 = %f +/- %f',coeffs(3),se(3))
