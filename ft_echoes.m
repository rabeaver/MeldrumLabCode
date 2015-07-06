%%
clear
clc
close all

real_data = load('Re_echoData.dat');
imag_data = load('Im_echoData.dat');
parfilestem = 'acqu';

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

T = params.acqTime/params.nrPts/1000;    %sample time (s)
Fs = 1/T;                           %sampling freq
L = params.nrPts;                   %length of signal (nrPts)
t = (1:L)*T;

%%
plot(t*1000,real_data(2,:))
xlabel('time (milliseconds)')

%%
NFFT = 2^nextpow2(L);
Y = fft(real_data(2,:),NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
plot(f/851.52,2*abs(Y(1:NFFT/2+1)));
xlabel('distance (um)')
