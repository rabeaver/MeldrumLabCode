clear
close all
clc


spectrometer = 'Kea'; %'Tecmag'
datadir = '/Volumes/ISC1026/Data/TKM/GCI_Acrylate/17May2016/Green_ManualProfile_FT/';
filenum = '1';
datafile = '/data';


pw = 2e-6;                      % hard pulse length
nPts = 16;                          % # of acqu points
omitPtsBack = 0;                    % the number of points at the end of each echo window that are zeros from the spectrometer
omitPtsFront = 0;                    % the number of points at the beginning of each echo window to zero
nEchoes = 256;                      % Echoes
omitEchoes = 0;                     % numner of echoes to remove from data
tD = 2e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 62;                           % us
noisePoints = 4;                    % number of points for measuring noise
cutRefPts = 0;                     %if necessary, can cut the data from the reference scan by half this value on each end of the acq window
                                   %use only if nPts for CHIRP on and CHIRP off expts don't match

zf = 2;                             % levels of zero filling
apodize = 0;                        % Gaussian apodization on (1) or off (0)?
apofac = 5;                         % Amount of Apodization


% ===================================
% === END User-defined paramaters ===
% ===================================

if strcmp(spectrometer,'Tecmag')==1;
    G = 6.59;                           % T m-1, B0 field gradient
elseif strcmp(spectrometer,'Kea')==1;    
    G = 23.87;
end

gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1

T = tD;                             % Sample time
Fs = 1/T;                           % Sampling frequency 
L = (nPts-omitPtsFront-omitPtsBack)*(2^zf);          % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

echoVec = tE*(omitEchoes+1):tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/(gamma*G);                    % um, 280.47 Hz/um (for PM25)

%% Import CHIRP data

for ii = 1:9
    if strcmp(spectrometer,'Tecmag')==1;
        [ap , spec] = readTecmag4d(strcat(datadir,num2str(ii),datafile,'.tnt'));
    elseif strcmp(spectrometer,'Kea')==1;
        [ap , spec] = readKea4d(strcat(datadir,num2str(ii),datafile,'.2d'));
    end
    
    dat = reshape(spec, nPts, nEchoes);
    dat = dat(1+omitPtsFront:end-omitPtsBack,omitEchoes+1:end);
    
    % Apodization, zero filling, do FFT
    pVec = 1:1:(nPts-omitPtsBack-omitPtsFront);
    filt = exp(-(pVec-(nPts-omitPtsBack-omitPtsFront)/2).^2/((nPts-omitPtsBack-omitPtsFront)/apofac)^2);
    filt = repmat(filt',1,nEchoes-omitEchoes);
    
    if apodize == 1
        dat = dat .* filt;
    end
    
    dat = padarray(dat, size(dat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's
    
    profile(:,:,ii) = flipud(fftshift(fft(dat,NFFT)/L, 1)); % Performs FFT algorithm
    z(ii,:) = f/(gamma*G)-(ii-1)*(100);                    % um, 280.47 Hz/um (for PM25)
end
    z = z';

%%    
figure(1)
for jj = 1:9
    subplot(3,3,jj)
    pcolor(echoVec,z(:,jj),abs(profile(:,:,jj))); shading flat
    caxis([0 max(max(max(abs(profile))))])
    ylim([-1050 250])
    xlim([0 nEchoes*tE])
end


