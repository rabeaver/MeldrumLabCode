clear
clc
close all

%%
filename = 'GdH2O_15mM_GreenVial_T1Sat_13July2015.tnt';
filedir = '/Users/jaredking/Documents/Chemistry/Research/CHIRP/';

fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tEcho = 500; %us
nEchoes = 16;
nPts = 67;
nPtsBlank = 3;
nT1Pts = 38;
T1min = 0.05; %ms
T1max = 5.437; %ms
T1guess = 1; %ms 

T1vector = linspace(T1min,T1max,nT1Pts); % Linspace T1sat
echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

% T1vector = logspace(log10(T1min),log10(T1max),nT1Pts); % Logspace T1sat

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


