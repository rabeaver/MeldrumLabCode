clear
clc
close all

data = load('/Users/tyler/Dropbox/Data/Biosensors/NPNa_BSA_Titration/BSA_Blank/T2_D/4/dataRe.dat');
echoVector = (60:60:1024*60)*1e-3;
deltaVec = logspace(log10(15e-6),log10(815e-6),7);
gamma = 267.513e6; %s-1 T-1
GMHz = 1016; %MHz/m
gammaMHz = 42.576; %MHz T-1
G = GMHz/gammaMHz; %T m-1
DELTA = 1e-3;
xaxis = -gamma^2.*G^2.*deltaVec.^2.*(DELTA+(2/3).*deltaVec);


dataD = sum(data(:,1:32),2);
dataD = dataD./dataD(1);
% [p,S,mu] = polyfit(xaxis',log(dataD),1)
cftool(xaxis(3:7),log(dataD(3:7)))
%%
% plot(data(5,:))

guesses = [0.4, 30, 0.6, 3];
CI = 90;

% [xfit,ypred,beta,beta_err,resid] = bidecay_t2fit(echoVector,data(5,:),guesses,CI);
% [xfit,ypred,beta,pm,resid] = monodecay_t2fit_simple(echoVector,data(1,:),guesses,CI);

figure(1)
subplot(2,1,1)
hold on
scatter(echoVector,data(1,:))
plot(xfit,ypred,'-r')
subplot(2,1,2)
scatter(echoVector,resid)

rmsd = sqrt(mean(resid.^2))