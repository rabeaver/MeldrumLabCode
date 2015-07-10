clear
clc
close all

%%
filename = 'DoubleSample_50mMand5mM_GdH2O_CPMG_10July2015.tnt';
filedir = 'C:\Users\NMRLab\Desktop\CHIRP\';

fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tEcho = 500; %us
nEchoes = 128;
nPts = 69;
nPtsBlank = 5;

echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

data = reshape(spec,nPts,nEchoes);
data = data(1:(nPts-nPtsBlank),:);
dataInt = sum(data,1);
dataIntRe = real(dataInt)./max(real(dataInt));
dataIntIm = imag(dataInt)./max(real(dataInt));

guess = [1 15e-03];
beta = nlinfit(echoVector,dataIntRe, @t2monofit_simple, guess);
ypred = t2monofit_simple(beta,echoVector);

figure(1)
hold on
plot(echoVector,dataIntRe);
plot(echoVector,dataIntIm);
plot(echoVector,ypred,'-r');