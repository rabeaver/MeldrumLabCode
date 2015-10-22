clear
clc
close all

%%

% Input filename, - .tnt
filename = '150mMGdH2O_T1IR_BURP__tE400u_10Sep2015';
filedir = '/Users/jaredking/Documents/Chemistry/Research/Raw Data/150mMGd/';

fileloc = strcat(filedir,filename,'.tnt');

% Read file
[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);

% Input experiment parameters
tEcho = 400; %us

nEchoes = 8;
nPts = 40;
nPtsBlank = 2;
nT1Pts = 21;
T1min = 0.1; %ms
T1max = 1.2; %ms
echoVector = (tEcho:tEcho:nEchoes*tEcho); % T2 vector


% Specify lin or log spaced points
linORlog = 0; % 0 for linearly space and 1 for log spaced

% Make T1vector
if linORlog == 0
    T1vector = linspace((T1min),(T1max),nT1Pts); % Linspace T1 points
else
    T1vector = logspace(log10(T1min),log10(T1max),nT1Pts); % Logspace T1sat
end

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
