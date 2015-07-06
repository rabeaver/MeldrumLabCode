%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is for processing data collected with the Tecmag using the
% dipolar-encoded longitudinal magnetization decay sequence
% (LMDecay_26Jan2015.tps). You must input:
% 1. directory of data (be sure to phase first)
% 2. the file name
% 3. number of echoes
% 4. echo time in us
% 5. the minimum and maximum evolution time (given from the 2d array in TNMR, 2*(t+pw))
% 6. the number of acuired points to use
% This will calculate a log-spaced vector of evolution points, reshape and 
% integrate the data, and plot the decaying signal intensity against both tau and tau^2.
% 
% TKM, 3 Feb 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all 
clc

%user input parameters
pathdir = 'C:\Users\maoplinger\Desktop\';
filename = 'Nitrile_LMDecay_03_10Feb';
% filename = 'Nitrile_LMDecay_03_10Feb';
nEchoes = 32;
tEcho  = 150; %us
tEvoMin = 500; %us
tEvoMax = 10000; %us
echoVector = tEcho*1e-3*(1:1:nEchoes);
nPtsUse = 64;

%load data, determine 2D points
cd(pathdir)
[ap,spec,spec2,~,~] = readTecmag4d(filename);
n2D = size(spec2,1);

%build tau axis (evolution time)
tEvo = logspace(log10(tEvoMin),log10(tEvoMax),n2D)*1e-3; %ms

%import, reshape, and sum data
nPts = size(spec2,2)/nEchoes;
data = reshape(spec2,n2D,nPts,nEchoes);
data = data(:,1:nPtsUse,:);
dataSum = sum(data,2);
dataSum = reshape(dataSum,n2D,nEchoes);
dataInt = sum(dataSum,2);

%plotting
figure(1)
surf(tEvo,echoVector,real(dataSum')); shading flat
xlabel('evolution time (ms)')
ylabel('echo vector (ms)')
 
figure(2)
subplot(1,2,1)
scatter(tEvo,real(dataInt)./max(real(dataInt)))
xlabel('\tau (ms)')
ylabel('normalized LM intensity (arb)')
ylim([0 1])
subplot(1,2,2)
scatter(tEvo.^2,real(dataInt)./max(real(dataInt)))
% xlim([0 10])
ylim([0 1])
xlabel('\tau^2 (ms^2)')
ylabel('normalized LM intensity (arb)')