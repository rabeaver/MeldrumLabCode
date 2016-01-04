clear
clc
close all

%% Get data

% Get Parameters
parfilestem = sprintf('%s', 'C:\Users\fjmorin\Desktop\PEGDA100%\100%PEGDA_180scure_slide_01Dec\3\acqu');

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
file = sprintf('%s', 'C:\Users\fjmorin\Desktop\PEGDA100%\100%PEGDA_180scure_slide_01Dec\3\data.csv');
data = load(file); % Open datafile

% Separate data
dataRe = data(1:params.nrEchoes,:);
% dataIm = data((params.nrEchoes+1):(params.nrEchoes*2),:);

% Time vector
echoVec = (1:params.nrEchoes)*params.echoTime/1e6;

%% User Defined Parameters

% Range tested in ILT
T2min = 6e-05; % Minimum T2 time (s)
T2max = 1e0; % Max T2 time (s)
% --

smoothing = 1e08; % Smoothing parameter "Alpha" (1e7 is a good place to start)
beta = -1; % Not really user defined... Never changes
steps = 256; % Number of points in ILT; can't be bigger than number of echoes

%% Do ILT

[spec, tau, chisq] = upnnlsmooth1D(dataRe(:,2), echoVec', T2min, T2max, smoothing, beta, steps);

%% Plot Spectrum

figure(1)
hold on
plot(tau, spec)
set(gca, 'xscale', 'log')
hold off
