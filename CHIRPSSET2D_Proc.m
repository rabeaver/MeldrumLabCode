clear
clc
close all

%
% CHIRP params
% ===================================
% ===== User-defined paramaters =====
% ===================================

Pchirp = 0.000397; % CHIRP Pulse Length (s)
pw     = 6e-6; %hard pulse length
sliceheight = 0.350; %mm

nPts = 42; % # of acqu points
omitPts = 0; %the number of points that are zeros from the spectrometer
nEchoes = 64; % Echoes
omitEchoes = 2; %numner of echoes to remove from data
tD = 5e-6; % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 300; %us

zf = 1;                             % levels of zero filling
apodize = 0;                        %Gaussian apodization on (1) or off (0)?
apofac = 5;                         % Amount of Apodization

deltaMax = 0.8e-3; % lil deltamax time in s
% deltaMin = deltaMax-2*Pchirp-pw/2; % calculates the minimum value of delta from the chirp
DELTA = 0.5e-3; % Big delta time in s

% ===================================
% === END User-defined paramaters ===
% ===================================

G = 6.59;                           % T m-1, B0 field gradient
gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % Hz T-1 radians
BWchirp = sliceheight*G*gamma*1000; % CHIRP bandwidth (Hz)

T = tD;                             % Sample time
Fs = 1/T;                           % Sampling frequency 
L = (nPts-omitPts)*(2^zf);      % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

echoVec = tE:tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/280.47;                       % um, 280.47 Hz/um (for PM25)

%

datadir = 'C:\CommonData\CHIRP\T2D\';
datafile = 'CHIRP_GdWaterSieves_T2D_800usd_500usD_397usCHIRP_350um_32pwr_4096sc_31Oct2015_result';


% Import CHIRP data
[~ , spec, spec2, ~] = readTecmag4d(strcat(datadir,datafile,'.tnt'));

% CHIRPdat = spec(1,:);
% spec = spec2(nnn, :);
CHIRPdat = reshape(spec, nPts, nEchoes);
CHIRPdat = CHIRPdat(1:end-omitPts,omitEchoes+1:end);

%

pVec = 1:1:(nPts-omitPts);
filt = exp(-(pVec-(nPts-omitPts)/2).^2/((nPts-omitPts)/apofac)^2);
filt = repmat(filt',1,nEchoes-omitEchoes);

if apodize == 1
    CHIRPdat = CHIRPdat .* filt;
end

CHIRPdat = padarray(CHIRPdat, size(CHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

T2Dprofiles = fftshift(fft(CHIRPdat,NFFT)/L, 1); % Performs FFT algorithm

figure(1)
subplot(1,2,1)
hold on
plot(t*1e6,real(CHIRPdat(:,1)));
xlabel('time [us]')

subplot(1,2,2)
hold on
plot(z,2*abs(T2Dprofiles(:,1)),'LineWidth',1.5);
xlabel('real space [um]')
title('Plot of first T1T2 FFT Profile and Echo')

hold off
figure(2)
hold on
surf(abs(T2Dprofiles));
shading flat;
title('Surface plot of T1T2 FFT Profiles')
hold off


% No CHIRP load section
% close all


noCHIRPfile = 'noCHIRP_GdWaterSieves_T2D_800usd_500usD_4096sc_31Oct2015_result';
[~,spec,spec2] = readTecmag4d(strcat(datadir,noCHIRPfile,'.tnt'));
data = reshape(spec,nPts,nEchoes);

% No CHIRP raw data and fft profiles
% data = spec2(2,:);
noCHIRPdat = reshape(data, nPts, nEchoes);
noCHIRPdat = noCHIRPdat(1:end-omitPts,omitEchoes+1:end);
if apodize == 1
    noCHIRPdat = noCHIRPdat .* filt;
end

noCHIRPdat = padarray(noCHIRPdat, size(noCHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's


CPprofiles = fftshift(fft(noCHIRPdat,NFFT)/L,1);

% Plot first T1-T2 profile and coil profile
% 
figure(3)
subplot(1,2,1)
hold on
plot(t*1e6,real(noCHIRPdat(:,1)));
xlabel('time [us]')
subplot(1,2,2)
hold on
plot(z,2*abs(CPprofiles(:,1)),'LineWidth',1.5);
xlabel('real space [um]')
title('Plot of first T1T2 FFT Profile and Echo')
hold off

figure(4)
hold on
surf(abs(CPprofiles));
shading flat;
title('Surface plot of T1T2 FFT Profiles')
hold off


% Plot first T1T2 profile and coil profile

% close all

figure(5)
hold on
plot(z,abs(CPprofiles(:,1))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','k')
plot(z,abs(T2Dprofiles(:,1))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','r')
hold off

xlabel('{\it z} (um)','fontsize',12)
title('T1-T2 & coil profiles')
set(gca,'Fontsize',12,'linewidth',2)

% Coil Sensitivity Correction

for k = 1:nEchoes-omitEchoes
    pcorr(:,k) = abs(CPprofiles(:,1));
end

T2Dprofcorr = T2Dprofiles./pcorr;

% figure(6)
% pcolor(abs(T1T2profcorr)); 
% colormap('jet');
% shading interp;
% colorbar('linewidth',2)
% caxis([0 1])
% title('Coil sensitivity corrected T1-T2 profiles')

% Find Optimal data range with these figures
% close all

figure(8)
plot(abs(T2Dprofiles(:,1)))

% t1_fig7 = linspace(deltaMax,0,NFFT);
t1_fig7=Pchirp*(BWchirp/2-f)/BWchirp;

%%
close all
figure(7)
subplot(2,1,1)
plot(abs(T2Dprofcorr(:,1:8)))
xlim([0 NFFT])
ylim([0 2])
subplot(2,1,2)
plot(t1_fig7,abs(T2Dprofcorr(:,1:8)))
line([0 0],[-2 2])
line([Pchirp Pchirp],[-2 2])
xlim([min(t1_fig7), max(t1_fig7)]);
ylim([0 2])
set(gca,'XDir','reverse')
xlabel('CHIRPtime (s)')



%% Data Range and Inversion


minind= 58;
maxind = 84;
% this is where I'm starting to put in some diffusion code. 


% tAxist2d = diff(t1_fig7); % gives the distance between each time point based on the time figure above
% tDift2d  = abs(round((deltaMax-deltaMin)/tAxist2d(1))); % makes into an integer based on the time domainabove

T2Ddat = abs(T2Dprofcorr(minind:maxind,:)); %crops data set according to above indices
deltaSteps = 2*t1_fig7(minind:maxind); %added a factor of 2 to account for how the refocusing makes the actual little-delta diffusion time range from 0 to the full delta time. (30 Oct 2015 TKM)

yD = log(T2Ddat./T2Ddat(1))';
xD = -gammaRad^2*G^2.*deltaSteps.^2.*(DELTA + (2/3)*deltaSteps);

T2Dsize = size(T2Ddat,1); % cuts down delta points to math those selected for the indices,assuming that the 
%first point is the first part of the chirp

figure(9)
plot(xD(1:T2Dsize),(yD(1,:))) 


% time axis are set depending on evenly spaced time points and assuming
% that the FIRST POINT taken for the indices is the 1st point of the chirp
% pulse

%% CF tool
a = 1;
b = T2Dsize;

cftool(xD(a:b),yD(1,a:b))


%% surf of all T1-T2 Profiles

figure(10)
surf(echoVec(omitEchoes+1:end)/1000,deltaSteps*1e6,T2Ddat);
shading flat
% xlabel('echo time [ms]')
% ylabel('delta [us]')
colormap('jet');
% shading interp;
colorbar 
xlabel('{\it T}_2 (ms)'); 
ylabel('{\it delta} (us)');
title('D-T2 data')


%% Save data, display ILT Data params
close all

T2Ddat2 = flipud(T2Ddat);
save(strcat(datadir,datafile, '.dat'), 'T2Ddat2', '-ascii')
size(T2Ddat)
1e6*abs(deltaSteps(1)-deltaSteps(end))
1e6*[min(deltaSteps), max(deltaSteps)]
