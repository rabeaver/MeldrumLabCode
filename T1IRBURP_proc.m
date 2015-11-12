clear
clc
close all

%%

% Input filename, - .tnt
filename = '15mMGdH2O_BigSample_T1IRLong_11Nov2015_result';
filedir = '/Users/tyler/Desktop/CHIRP_Manuscript/Raw Data/11Nov2015_BigGdWater_512Scans/';
fileloc = strcat(filedir,filename,'.tnt');

% Read file
[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);

% Input experiment parameters

tEcho = 700; %us
nEchoes = 16;

echoVector = (tEcho:tEcho:nEchoes*tEcho); % T2 vector


nPts = 76;
nPtsBlank = 4;
nT1Pts = 21;
T1min = 0.05; %ms
T1max = 24.95; %ms
noisePoints = 10; %number of points to use for noise at beginning and end of each acqu period


% Specify lin or log spaced points
linORlog = 0; % 0 for linearly space and 1 for log spaced

% Make T1vector
if linORlog == 0
    T1vector = linspace((T1min),(T1max),nT1Pts); % Linspace T1 points
else
    T1vector = logspace(log10(T1min),log10(T1max),nT1Pts); % Logspace T1sat
end

%% SNR calc
data = reshape(spec2,nT1Pts,nPts,nEchoes);
data = data(:,1:nPts-nPtsBlank,:);
n1 = data(:,1:noisePoints,:);
n2 = data(:,nPts-nPtsBlank-noisePoints:end,:);
ndata = cat(2,n1,n2);
ndata = reshape(ndata,nT1Pts,(2*noisePoints+1)*nEchoes);
sdata = reshape(data,nT1Pts,(nPts-nPtsBlank)*nEchoes);

% data = reshape(data,nT1Pts,(nPts-nPtsBlank)*nEchoes);

%
S = max(abs(sdata(end,:)));
N = rms(ndata(end,:));

% SNR = snr(S,N)
SNR = S/N

% figure
% hold on
% plot(S)
% plot(N)
%% Make 2D data set for T1IRT2 ILT

data = reshape(spec2',nPts,nEchoes,nT1Pts);
data = data(1:(nPts-nPtsBlank),:,:);
data2d = sum(real(data),1);
data2d = reshape(data2d,nEchoes, nT1Pts);
data2d = data2d';

% Plot of data
surf(data2d); shading flat

% Save data in specified directory with the same filename and ".dat"
% extension
save(strcat(filedir,filename,'.dat'), 'data2d', '-ascii')

%% 1D Fits

%T1 (A*(1-2*exp(-x/T1))
cftool(T1vector, data2d(:,1)./max(data2d(:,1)));

%T2 (A*exp(-x/T2))
cftool(echoVector/10^3, data2d(end,:)'./max(data2d(end,:)))
