clear
clc
close all

%%
filename = '1500uM_GdH2O_BeadPack_CPMG_13July2015.tnt';
filedir = 'C:\Users\NMRLab\Desktop\CHIRP\';

fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tEcho = 500; %us
nEchoes = 32;
nPts = 69;
nPtsBlank = 5;

echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

data = reshape(spec,nPts,nEchoes);
data = data(1:(nPts-nPtsBlank),:);
dataInt = sum(data,1);
dataIntRe = real(dataInt)./max(real(dataInt));
dataIntIm = imag(dataInt)./max(real(dataInt));

guess = [0.3 2e-03 0.6 6e-03];
beta = nlinfit(echoVector,dataIntRe, @t2bifit_simple, guess);
ypred = t2bifit_simple(beta,echoVector);

figure(1)
hold on
plot(echoVector,dataIntRe);
plot(echoVector,dataIntIm);
plot(echoVector,ypred,'-r');