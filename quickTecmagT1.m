clear
clc
close all

%%

filename = 'Membrane_w_pureH2O_T1IR_BURP_28June2016.tnt';
filedir = 'C:\CommonData\Membranes\PureWater\';

fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tEcho = 150; %us

nEchoes = 64;
nPts = 16;
nPtsBlank = 0;
nT1Pts = 11;
T1min = 0.1; %ms
T1max = 200; %ms
T1guess = 50; %ms 

T1vector = linspace((T1min),(T1max),nT1Pts); % Linspace T1sat
echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

% T1vector = logspace(log10(T1min),log10(T1max),nT1Pts); % Logspace T1sat

data = reshape(spec2',nPts,nEchoes,nT1Pts);
data = data(1:(nPts-nPtsBlank),:,:);
dataInt = sum(sum(data,1),2);
dataInt = reshape(dataInt,1,nT1Pts);
dataInt = dataInt./(dataInt(1));
dataIntRe = real(dataInt);
dataIntIm = imag(dataInt);


%% cftool
cftool(T1vector, -dataIntRe)

%%
% guesses = [1, max(dataIntRe), T1guess];
[beta, Resids, J] = nlinfit(T1vector,dataIntRe,@inversionRecoverySimple,T1guess);
ypred = inversionRecoverySimple(beta,T1vector);
CI = nlparci(beta, Resids, 'jacobian', J);


figure(1)
hold on
scatter(T1vector,dataIntRe);
plot(T1vector,ypred);
