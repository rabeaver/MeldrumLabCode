clear
clc
close all

%%
% CHIRP params
% ===================================
% ===== User-defined paramaters =====
% ===================================

Pchirp = 300e-06; % CHIRP Pulse Length (s)

sliceheight = 0.050; %mm
PreCPMGdelay = 40e-6; %s


nPts = 128; % # of acqu points
nEchoes = 16; % Echoes
tD = 3e-6; % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 500; %us
nScans = 32; %number of scans

omitEchoes = 0; %the number of echoes to skip
noisePoints = 3; %number of points at beginning and end of each acqu period for noise
omitPts = 0; %blank spectrometer points to skip
zf = 1;                             % levels of zero filling

%pwrRange = [80 70 60 50 40 30 20 10 6]; % absolute value of range of CHIRP powers being evaluated
pwrRange = [40 39 38 37 36 35 34 33 32 31 30];
%pwrRange = [39 38 31];
SNR = [];
SNR_perRtScans = [];
maxsig = [];
% ===================================
% === END User-defined paramaters ===
% ===================================

G = 23.87;
    
gamma = 42.576;                     % MHz T-1
BWchirp = sliceheight*G*gamma*1000; % CHIRP bandwidth (Hz)

T = tD;                             % Sample time
Fs = 1/T;                           % Sampling frequency 
L = nPts*(2^zf);                    % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

echoVec = tE:tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/(gamma*G);                    % um, 280.47 Hz/um (for PM25)

%%
% Import data
for i = 1:length(pwrRange)
    
spectrometer = 'Kea'; %'Tecmag' or 'Kea'
datadir = 'C:\Users\jnking01\Desktop\50um_sliceheight\400us_cpw\40to30dBRange\';
datafile = strcat('400us_neg', num2str(pwrRange(i)), 'db\1\data');
noCHIRPfile = strcat('400us_neg', num2str(pwrRange(i)), 'db\1\data');
filenameExt = '';

if strcmp(spectrometer,'Tecmag')==1;
    [ap , spec] = readTecmag4d(strcat(datadir,datafile,'.tnt'));
elseif strcmp(spectrometer,'Kea')==1;
    [ap , spec] = readKea4d(strcat(datadir,datafile,'.2d'));
end

% CHIRPdat = spec(1,:);
% spec = spec2(nnn, :);
CHIRPdat = reshape(spec, nPts, nEchoes);
%CHIRPdat = CHIRPdat(1:(end-omitPts),omitEchoes+1:end);

%%

%% Find Maximum
maxsig(i) = max(max(abs(CHIRPdat)));

%% SNR calc (two sections)
n1 = CHIRPdat(1:noisePoints,:);
n2 = CHIRPdat(nPts-noisePoints:end,:);
n = cat(1,n1,n2);
n = reshape(n,1,(2*noisePoints+1)*(nEchoes-omitEchoes));
s = reshape(CHIRPdat,1,(nPts-omitPts)*(nEchoes-omitEchoes));

S = max(abs(s));
N = rms(n);

SNR(i) = S/N;
SNR_perRtScans(i) = SNR(i)/sqrt(nScans); % Set up for power cal, aka only CHIRP on, no CHIRP-less scan

%%
% figure()
% hold on
% plot(echoVec,sum(abs(CHIRPdat),1));
% title(strcat('CPMG at power -',num2str(pwrRange(i))))
% hold off

realDat = real(CHIRPdat);
imagDat = imag(CHIRPdat);

figure()
hold on
plot(realDat(:,2))
plot(imagDat(:,2))
ylim([-1.5 1.5])
title(strcat('Echo 2 at power -', num2str(pwrRange(i))))
hold off

end

% figure(1)
% plot(-pwrRange,maxsig, 'r')
% title('Maximum signal vs. power')
% 
% 
% % figure(2)
% % plot(-pwrRange,SNR_perRtScans, 'g')
% % title('SNR per root scan vs. power')
% 
% figure(3)
% plot(-pwrRange,SNR,'b')
% title('SNR vs. power')