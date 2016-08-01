clear
clc
close all

% USER-DEFINED PARAMETERS
filename = 'data.2d';
filedir = 'Z:\Data\VJL\Coatings\CPMG\PEGDA_100%_11s_10%_30June2016\1\';

omitEchoes = 2; %front-end echoes to omit
G = 23.87; %T/m
zf = 1; %levels of zeo-filling
% END USER-DEFINED PARAMETERS


fileloc = strcat(filedir,filename);
parloc  = strcat(filedir,'acqu.par');

[ap,spec] = readKea4d(fileloc);
tE = readpar_Kea(parloc,'echoTime')*1e-6;
tD = readpar_Kea(parloc,'dwellTime')*1e-6;
nrEchoes = ap.yDim;
nrPts = ap.xDim;

echoVec = (omitEchoes+1)*tE:tE:nrEchoes*tE;
spec2d = reshape(spec,nrPts,nrEchoes);
spec2d = spec2d(:,omitEchoes+1:end);

gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1

T = tD;                             % Sample time
Fs = 1/T;                           % Sampling frequency 
L = (nrPts)*(2^zf);          % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

t = (-(L-1)/2:L/2)*T;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/(gamma*G);                    % um, 280.47 Hz/um (for PM25)


FTdat = padarray(spec2d, size(spec2d(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

FTT2 = (fftshift(fft(FTdat,NFFT)/L, 1)); % Performs FFT algorithm

% Plot CHIRP results
figure(1)
subplot(1,2,1)
plot(t*1e6,real(FTT2(:,1)));
xlabel('time [us]')
subplot(1,2,2)
plot(z,2*abs(FTT2(:,1)),'LineWidth',1.5);
xlabel('real space [um]')
title('Plot of first T2-D FFT Profile and Echo')

figure(2)
surf(echoVec'*1e3,z,abs(FTT2));
shading flat;
title('Surface plot of T2-D FFT Profiles')
xlabel('T2 [ms]')
ylabel('z [um]')
view([0 90])

[C,I] = max(abs(FTT2));

figure(3)
scatter(echoVec,abs(FTT2(I(1),:)));
% guess = [0.4 1 7]; %T2 in ms
% [beta,R,J,CovB] = nlinfit(1e-3*spec(:,1),spec(:,2)./spec(1,2), @t2bifit_ampSumFixed, guess);
% ci = nlparci(beta,R,'jacobian',J);
% 
% ypred = t2bifit_ampSumFixed(beta,1e-3*spec(:,1));
% 
% hh = figure(1);
% hold on
% scatter(1e-3*spec(:,1),spec(:,2)./spec(1,2));
% plot(1e-3*spec(:,1),ypred,'-r');
% xlabel('time/ms')
% % pubgraph(hh, 16, 2, 'w')
% 
% sprintf('T2 = %f +/- %.1g ms.',beta(2), beta(2)-ci(2,1))