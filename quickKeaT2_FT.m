clear
clc
close all

% USER-DEFINED PARAMETERS
filename = 'data.2d';
filedir = 'C:\CommonData\JNK\UFT2D\PaintForNick\M212_CPMG_23June2016\2\';

zf = 1;

omitEchoes = 0; %front-end echoes to omit
% END USER-DEFINED PARAMETERS

G = 23.87;
gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1

fileloc = strcat(filedir,filename);
parloc  = strcat(filedir,'acqu.par');

[ap,spec] = readKea4d(fileloc);
tE = readpar_Kea(parloc,'echoTime')*1e-6;
nrEchoes = ap.yDim;
nrPts = ap.xDim;
acqTime = readpar_Kea(parloc,'acqTime');

echoVec = (omitEchoes+1)*tE:tE:nrEchoes*tE;
spec2d = reshape(spec,nrPts,nrEchoes);
spec2d = spec2d(:,omitEchoes+1:end);

T = acqTime/nrPts;                             % Sample time
Fs = 1/T;                           % Sampling frequency 
L = (nrPts)*(2^zf);          % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

t = (-(L-1)/2:L/2)*T;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/(gamma*G);                    % um, 280.47 Hz/um (for PM25)

dat = padarray(spec2d, size(spec2d(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

profiles = flipud(fftshift(fft(dat,NFFT)/L, 1)); % Performs FFT algorithm

plot(z,abs(profiles(:,1)));