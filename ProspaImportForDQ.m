clear
clc
close all

%% Get data

% Get Parameters
filedir = sprintf('%s', 'C:\Users\jhyu\Desktop\1\');
parfilestem = strcat(filedir,'acqu');

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
dataRe = load(strcat(filedir,'dataRe.dat')); % Open datafile
dataIm = load(strcat(filedir,'dataIm.dat'));

% Time vector
tau = load(strcat(filedir,'data.dat'));
tau = tau(:,1);
echoVec = (1:params.nrEchoes)*params.echoTime;

dataRe2 = reshape(dataRe,length(tau),params.nrPts,params.nrEchoes);
dataRe2 = sum(dataRe2,2);
dataRe2 = reshape(dataRe2,length(tau),params.nrEchoes);

dataIm2 = reshape(dataIm,length(tau),params.nrPts,params.nrEchoes);
dataIm2 = sum(dataIm2,2);
dataIm2 = reshape(dataIm2,length(tau),params.nrEchoes);

dataCp = complex(dataRe2,dataIm2);

%%
plot(echoVec,dataRe2(1,:))

%%
surf(echoVec,tau,real(dataCp)); shading flat