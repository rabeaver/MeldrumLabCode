%% BIEXP FIT
clear 
close all
clc

datadir = '/Users/tyler/Desktop/Bead Pack/1500uMGd_BEADPACK_50u/1/';
datafile = 'data2.csv';
parfilestem = strcat(datadir,'acqu');

% Load data and params
data = load(strcat(datadir,datafile));
params.acqTime = readpar_Kea(strcat(parfilestem,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(parfilestem,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(parfilestem,'.par'),'nrScans');
params.rxPhase = readpar_Kea(strcat(parfilestem,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(parfilestem,'.par'),'rxGain');
params.nrPts = readpar_Kea(strcat(parfilestem,'.par'),'nrPnts');
params.repTime = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.b1Freq = readpar_Kea(strcat(parfilestem,'.par'),'b1Freq');
params.nrEchoes = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoes');
params.echoTime = readpar_Kea(strcat(parfilestem,'.par'),'echoTime');

echoVector = params.echoTime:params.echoTime:params.nrEchoes*params.echoTime;
echoVector = echoVector'*1e-3; %switch to ms

% reshape data for formatting
data = reshape(data,(params.nrEchoes*2),params.nrPts);
dataRe = data(1:params.nrEchoes,:);
dataIm = data((params.nrEchoes+1):(params.nrEchoes*2),:);
dataCp = complex(dataRe,dataIm);
dataInt = sum(dataCp,2);
dataInt = dataInt./max(dataInt);

% plot(echoVector,real(dataInt))

% Fit to exp decay
A_guess = 0.3;
t2_guess = 3;
A_guess2 = 0.7;
t2_guess2 = 10; %ms

guesses = [A_guess;t2_guess;A_guess2;t2_guess2];

CI = 90; %desired confidence interval in percent

for i = 1:10
    [xfit,ypred,coeffs,coeffs_err,residuals] = bidecay_t2fit_simple(echoVector,abs(dataInt),guesses,CI);
    guesses = coeffs;
end

figure
subplot(2,1,1)
hold on
plot(echoVector,abs(dataInt))
plot(xfit,ypred)
subplot(2,1,2)
plot(echoVector,residuals)

%% MONOEXP FIT
clear 
close all
clc

datadir = '/Users/tyler/Desktop/Bead Pack/1500uMGd_SOLUTION_80u/1/';
datafile = 'data2.csv';
parfilestem = strcat(datadir,'acqu');

% Load data and params
data = load(strcat(datadir,datafile));
params.acqTime = readpar_Kea(strcat(parfilestem,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(parfilestem,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(parfilestem,'.par'),'nrScans');
params.rxPhase = readpar_Kea(strcat(parfilestem,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(parfilestem,'.par'),'rxGain');
params.nrPts = readpar_Kea(strcat(parfilestem,'.par'),'nrPnts');
params.repTime = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.b1Freq = readpar_Kea(strcat(parfilestem,'.par'),'b1Freq');
params.nrEchoes = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoes');
params.echoTime = readpar_Kea(strcat(parfilestem,'.par'),'echoTime');

echoVector = params.echoTime:params.echoTime:params.nrEchoes*params.echoTime;
echoVector = echoVector'*1e-3; %switch to ms

% reshape data for formatting
data = reshape(data,(params.nrEchoes*2),params.nrPts);
dataRe = data(1:params.nrEchoes,:);
dataIm = data((params.nrEchoes+1):(params.nrEchoes*2),:);
dataCp = complex(dataRe,dataIm);
dataInt = sum(dataCp,2);
dataInt = dataInt./max(dataInt);

% plot(echoVector,real(dataInt))

% Fit to exp decay
A_guess = 0.5;
t2_guess = 5; %ms

guesses = [A_guess;t2_guess];

CI = 90; %desired confidence interval in percent
for i=1:10
    [xfit,ypred,coeffs,coeffs_err,residuals] = monodecay_t2fit_simple(echoVector,abs(dataInt),guesses,CI);
    guesses=coeffs;
end

figure
subplot(2,1,1)
hold on
plot(echoVector,abs(dataInt))
plot(xfit,ypred)
subplot(2,1,2)
plot(echoVector,residuals)
