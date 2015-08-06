clear
clc
close all

%% Get data

% Get Parameters
parfilestem = sprintf('%s', 'C:\Users\jhyu\Desktop\1\acqu');

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

% Datafile
file = sprintf('%s', 'C:\Users\jhyu\Desktop\1\dataRe.dat');
data = load(file); % Open datafile

% Separate data
dataRe = data(1:params.nrEchoes,:);
dataIm = data((params.nrEchoes+1):(params.nrEchoes*2),:);

% Time vector
echoVec = (1:params.nrEchoes)*params.echoTime;

%% User Defined Parameters

% Range tested in ILT
T2min = 1e-04; % Minimum T2 time (s)
T2max = 1e0; % Max T2 time (s)
% --

smoothing = 1e07; % Smoothing parameter "Alpha" (1e7 is a good place to start)
beta = -1; % Not really user defined... Never changes
steps = 64; % Number of points in ILT; can't be bigger than number of echoes

%% Do ILT

[spec, tau, chisq] = upnnlsmooth1D(dataRe', echoVec', T2min, T2max, smoothing, beta, steps);

%% Plot Spectrum

figure(1)
hold on
plot(spec, tau)
set('gca', 'xscale', 'log')
hold off
