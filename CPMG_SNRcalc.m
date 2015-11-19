clear
clc
close all

nScans = 1024;
%
Sfilename = strcat('GdWater_CPMG_',int2str(nScans),'scans_result.tnt');
% Nfilename = 'Glycerol_CPMG_16scans_12srep_NOISE_result.tnt';
filedir = '/Users/tyler/Dropbox/Data/SNRCheck/';
Sfileloc = strcat(filedir,Sfilename);
% Nfileloc = strcat(filedir,Nfilename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(Sfileloc);
% [ap,Nspec,spec2,spec3,spec4] = readTecmag4d(Nfileloc);
tEcho = 160; %us
nEchoes = 32;
nPts = 76;
nPtsBlank = 4;
nNoisePts = 10; %number of points at beginning and end of echo to measure noise

%
data = reshape(spec,nPts,nEchoes);
data = data(1:nPts-nPtsBlank,:);
Sspec = reshape(data,1,(nPts-nPtsBlank)*nEchoes);
N1 = data(1:nNoisePts,:);
N2 = data(nPts-nPtsBlank-nNoisePts:end,:);
Nspec = cat(1,N1,N2);
Nspec = reshape(Nspec,1,(2*nNoisePts+1)*nEchoes);
% a = real(mean(Sspec))
% b = var(Sspec)
% c = sqrt(b)
% 
% Sspec = sum(spec2,1);

SNR = max(abs(Sspec))/rms(Nspec)
%%

% echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;
% 
% data = reshape(spec,nPts,nEchoes);
% data = data(1:(nPts-nPtsBlank),:);
% dataInt = sum(data,1);
% dataIntRe = real(dataInt)./max(real(dataInt));
% dataIntIm = imag(dataInt)./max(real(dataInt));

% SNR calc
useEcho = 1;
Spoint = (useEcho+0.5)*nPts;
[~,Spoint] = max(abs(real(Sspec)));
S = (real(Sspec(Spoint-nPts/4:Spoint+nPts/4)));
Na = (imag(Sspec(Spoint-nPts/2:Spoint-nPts/4)));
Nb = (imag(Sspec(Spoint+nPts/4:Spoint+nPts/2-1)));

N = [Na, Nb];

% N = real(Sspec(Spoint-nPts/2:Spoint+nPts/2-4));
%
figure
hold on
% plot(S)
plot(N)

var(N)
mean(N)

%
SNR = snr(S,N)

%%
scans = [4 16 32 128 256 512 1024];
sscans = sqrt(scans);
SNR = [12.0325 18.5149 25.0474 33.1063 37.2905 40.3618 41.7243];

scatter(log2(scans),SNR)

%%
cftool(sscans,SNR)

