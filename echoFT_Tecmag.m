clear
clc
close all

%%
% ===================================
% ===== User-defined paramaters =====
% ===================================
%


datadir = 'C:\CommonData\TAMU\membrane holder_CuWater & Water test\';
datafile = 'Water_CPMG_13June2016_1m';


nPts = 30;                          % # of acqu points
nEchoes = 64;                      % Echoes% omitEchoes = 0;                     % numner of echoes to remove from data
tD = 6e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 450;                           % us

zf = 2;                             % levels of zero filling

% ===================================
% === END User-defined paramaters ===
% ===================================

G = 6.59;                           % T m-1, B0 field gradient
gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1

T = tD;                             % Sample time
Fs = 1/T;                           % Sampling frequency 
L = nPts*(2^zf);          % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

echoVec = tE:tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/(gamma*G);                       % um, 280.47 Hz/um (for PM25)

%% Import CHIRP data
[ap , spec] = readTecmag4d(strcat(datadir,datafile,'.tnt'));
dat = reshape(spec, nPts, nEchoes);

%% Zero filling, do FFT, plot
dat = padarray(dat, size(dat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

profile = (fftshift(fft(dat,NFFT)/L, 1)); % Performs FFT algorithm

figure(1)
surf(echoVec*1e-3,z,abs(profile)); shading flat
view([0 90])
xlabel('time [ms]')
ylabel('position [um]')