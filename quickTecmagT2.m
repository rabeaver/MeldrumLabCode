clear
clc
close all

%%
filename = 'Rubberstopper_T2_01_17Mar.tnt';
filedir = 'C:\Users\tkmeldrum\Desktop\';

fileloc = strcat(filedir,filename);

%%%%%%%%
[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tEcho = 150; %us
nEchoes = 128;
nPts = 69;
nPtsBlank = 5;
%%%%%%%%


echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

data = reshape(spec,nPts,nEchoes);
data = data(1:(nPts-nPtsBlank),:);
dataInt = sum(data,1);
dataIntRe = real(dataInt)./max(real(dataInt));
dataIntIm = imag(dataInt)./max(real(dataInt));

guess = nEchoes/10*tEcho*1e-6;
beta = nlinfit(echoVector,dataIntRe,@t2monofit_simple,guess);
ypred = t2monofit_simple(beta,echoVector);

figure(1)
hold on
scatter(echoVector,dataIntRe,'ok');
scatter(echoVector,dataIntIm,'or');
plot(echoVector,ypred,'-r');