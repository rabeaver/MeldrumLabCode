clear
clc
close all

%%
filename = '0P2_T1Sat_12.tnt';
filedir = 'C:\Users\jhyu\Desktop\';

fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tEcho = 150; %us
nEchoes = 512;
nPts = 35;
nPtsBlank = 2;
nT1Pts = 40;
T1min = .05; %ms
T1max = 1650; %ms
T1guess = 330; %ms 

T1vector = logspace(log10(T1min),log10(T1max),nT1Pts); % Linspace T1sat
echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

% T1vector = logspace(log10(T1min),log10(T1max),nT1Pts); % Logspace T1sat

data = reshape(spec2',nPts,nEchoes,nT1Pts);
data = data(1:(nPts-nPtsBlank),:,:);
dataInt = sum(sum(data,1),2);
dataInt = reshape(dataInt,1,nT1Pts);
dataIntRe = real(dataInt);
dataIntIm = imag(dataInt);

guesses = [1, max(dataIntRe), T1guess];
[beta, Resids, J] = nlinfit(T1vector,dataIntRe,@T1_recovery,guesses);
ypred = T1_recovery(beta,T1vector);
CI = nlparci(beta, Resids, 'jacobian', J);


figure(1)
hold on
scatter(T1vector,dataIntRe);
plot(T1vector,ypred);

%% Make 2D data set for T1SRT2

data2d = sum(abs(data),1);
data2d = reshape(data2d,nEchoes, nT1Pts);

surf(data2d); shading flat


