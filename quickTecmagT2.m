clear
clc
close all

%%
filename = 'CuWater_CPMG_21Apr2016_1.tnt';
filedir = 'C:\CommonData\TAMU\CuWater\';
fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);

tEcho = 400; %us
nEchoes = 64;
nPts = 54;
nPtsBlank = 4;

%% SNR calc


[~,Spoint] = max(real(spec2));
Spoint = Spoint + 3.5*nPts;
S = (real(spec(Spoint-nPts/2:Spoint+nPts/2)));
N = (imag(spec(Spoint-nPts/2:Spoint+nPts/2)));
% N = (real(specN(Spoint-nPts/2:Spoint+nPts/2)))';

SNR = snr(S,N)

figure
hold on
plot(S)
plot(N)

%%


echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

data = reshape(spec,nPts,nEchoes);
data = data(1:(nPts-nPtsBlank),:);
dataInt = sum(data,1);
dataIntRe = real(dataInt)./max(real(dataInt));
dataIntIm = imag(dataInt)./max(real(dataInt));

guess = [1 15e-3];% 0.6 6e-03];
beta = nlinfit(echoVector,dataIntRe, @t2monofit_simple, guess);
ypred = t2monofit_simple(beta,echoVector);

figure(1)
hold on
plot(echoVector,dataIntRe);
plot(echoVector,dataIntIm);
plot(echoVector,ypred,'-r');
xlabel('time/s')