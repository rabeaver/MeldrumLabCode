clear
clc
close all

filedir = '/Users/tyler/Desktop/Bead Pack/1500uMGd_BEADPACK_90u/1/';
filename = 'data.csv';
parname = 'acqu';
alpha = 1e10;
omitpoints = 0;

data = load(strcat(filedir,filename));

realData = data(omitpoints+1:end,2)./max(data(omitpoints+1:end,2));

params.acqTime = readpar_Kea(strcat(filedir,parname,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(filedir,parname,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(filedir,parname,'.par'),'nrScans');
params.rxPhase = readpar_Kea(strcat(filedir,parname,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(filedir,parname,'.par'),'rxGain');
params.nrPts = readpar_Kea(strcat(filedir,parname,'.par'),'nrPnts');
params.repTime = readpar_Kea(strcat(filedir,parname,'.par'),'repTime');
params.b1Freq = readpar_Kea(strcat(filedir,parname,'.par'),'b1Freq');
params.nrEchoes = readpar_Kea(strcat(filedir,parname,'.par'),'nrEchoes');
params.echoTime = readpar_Kea(strcat(filedir,parname,'.par'),'echoTime');

echoVector = params.echoTime:params.echoTime:params.echoTime*params.nrEchoes;
echoVector = echoVector*1e-6; %convert to s
lowLim = 10^-4; %
hiLim = 10^-1; %;
nrILTSteps = length(echoVector);

scatter(echoVector,realData)

[spectrum,tau,chisq,~     ] = upnnlsmooth1D(realData,echoVector',  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
semilogx(tau,spectrum)
xlabel('T_2 time [s]')
ylabel('intensity [arb]')
% ILTout = [tau',spectrum'];

% save('ILTout.dat','ILTout','-ascii');

