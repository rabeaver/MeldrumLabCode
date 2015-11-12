clear
clc
close all

%%

% Input filename, - .tnt
filename = '15mMGdH2O_T1IRLong_17Oct2015';
filedir = '/Users/tyler/Desktop/CHIRP_Manuscript/Raw Data/22Oct2015_1024scanProcessing/';
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
noisePts = 10; %number of points at front and back of each aqu period to use as noise
T1min = 0; %ms
T1max = 15000; %ms



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
N1 = data(:,1:noisePts,:);
N2 = data(:,nPts-noisePts-nPtsBlank:end,:);
dataN = cat(2,N1,N2);
noiseData = reshape(dataN,nT1Pts,(noisePts*2+1)*nEchoes);
%%
S = max(abs(spec2(1,:)));
N = rms(noiseData(1,:));


SNR = snr(S,N)
SNR2 = S/N


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
