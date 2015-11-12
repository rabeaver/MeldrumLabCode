clear
clc
close all

%
filename = 'Glycerol_CPMG_1024scans_12srep_result.tnt';
filedir = '/Users/tyler/Dropbox/Data/SNRCheck/';
fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tEcho = 700; %us
nEchoes = 128;
nPts = 76;
nPtsBlank = 0;

echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

data = reshape(spec,nPts,nEchoes);
data = data(1:(nPts-nPtsBlank),:);
dataInt = sum(data,1);
dataIntRe = real(dataInt)./max(real(dataInt));
dataIntIm = imag(dataInt)./max(real(dataInt));

% SNR calc
useEcho = 1;
Spoint = (useEcho+0.5)*nPts;
% [~,Spoint] = max(abs(real(spec)));
S = abs(real(spec(Spoint-nPts/2:Spoint+nPts/2)));
N = abs(imag(spec(Spoint-nPts/2:Spoint+nPts/2)));

figure
hold on
plot(S)
plot(N)

%
SNR = snr(S,N)
