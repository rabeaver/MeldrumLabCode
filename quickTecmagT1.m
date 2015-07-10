clear
clc
close all

%%
filename = 'GdH2O_50mMvial_T1Sat_10July2015.tnt';
filedir = 'C:\Users\NMRLab\Desktop\CHIRP\';

fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tEcho = 150; %us
nEchoes = 8;
nPts = 69;
nPtsBlank = 5;
nT1Pts = 11;
T1min = 0.05; %ms
T1max = 8; %ms
T1guess = 1.5; %ms 

T1vector = linspace(T1min,T1max,nT1Pts);
echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

T1vector = logspace(log10(T1min),log10(T1max),nT1Pts);

data = reshape(spec2',nPts,nEchoes,nT1Pts);
data = data(1:(nPts-nPtsBlank),:,:);
dataInt = sum(sum(data,1),2);
dataInt = reshape(dataInt,1,nT1Pts);
dataIntRe = real(dataInt);
dataIntIm = imag(dataInt);

guesses = [0, max(dataIntRe), T1guess];
beta = nlinfit(T1vector,dataIntRe,@T1_recovery,guesses);
ypred = T1_recovery(beta,T1vector);

figure(1)
hold on
scatter(T1vector,dataIntRe);
plot(T1vector,ypred);


