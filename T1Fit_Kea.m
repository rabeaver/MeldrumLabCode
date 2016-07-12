clear
clc
close all

%% T1 Fit from Kea Data

parfilestem = sprintf('%s', 'C:\Users\jnking01\Desktop\Atlas in ARLO 9.2\T1\1\acqu');

params.acqTime = readpar_Kea(strcat(parfilestem,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(parfilestem,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(parfilestem,'.par'),'nrScans');
params.rxPhase = readpar_Kea(strcat(parfilestem,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(parfilestem,'.par'),'rxGain');
params.nrPts = readpar_Kea(strcat(parfilestem,'.par'),'nrPnts');
params.repTime = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.repTime = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.b1Freq = readpar_Kea(strcat(parfilestem,'.par'),'b1Freq');
params.nrEchoes = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoes');
params.echoTime = readpar_Kea(strcat(parfilestem,'.par'),'echoTime');


file = sprintf('%s', 'C:\Users\jnking01\Desktop\Atlas in ARLO 9.2\T1\1\data.dat');
data = load(file);

data = data';
T1amp = data(1,:);
T1vector = data(2,:);


%% T1 Bifit

% Input guess, do fitting
guess = [0.02 1 100e-3];
[beta,Resids,J,covB] = nlinfit(T1amp,T1vector,@T1_recovery,guess);

% Confidence Interval and Margin of Error (+/-)
CI=nlparci(beta,Resids,'jacobian',J);
MOE_T1 = (CI(2,2) - CI(2,1))/2;

% Curve from fit
ypred = T1_recovery(beta, T1amp);

% plot fit
figure(1)
hold on
plot(T1amp,T1vector,'k*' )
plot(T1amp, ypred)


