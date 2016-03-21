clear
clc
close all

%%
% CHIRP params
% ===================================
% ===== User-defined paramaters =====
% ===================================
%


datadir = 'C:\CommonData\Acetone\';
datafile = 'AcetoneLarge_chirpSTE_short_18Mar2016_3_result';
noCHIRPfile = 'AcetoneLarge_nochirpSTE_short_18Mar2016_3_result';


Pchirp = 296.8e-6;                  % CHIRP Pulse Length (s)
pw     = 6e-6;                      % hard pulse length
sliceheight = 0.150;                % mm
rampPct = 0.01;                     % percent for the CHIRP power ramp to reach pMax

nPts = 54;                          % # of acqu points
omitPts = 4;                        % the number of points that are zeros from the spectrometer
nEchoes = 64;                      % Echoes
omitEchoes = 0;                     % numner of echoes to remove from data
tD = 6e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 400;                           % us
preCHIRPdelay = 0.2e-6;             % s
noisePoints = 2;                    % number of points for measuring noise

zf = 1;                             % levels of zero filling
apodize = 0;                        % Gaussian apodization on (1) or off (0)?
apofac = 5;                         % Amount of Apodization

delta = 0.6e-3;                       % little delta time (s)
DELTA = 0.5e-3;                       % Big delta time in s

% ===================================
% === END User-defined paramaters ===
% ===================================

G = 6.59;                           % T m-1, B0 field gradient
gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1
BWchirp = sliceheight*G*gamma*1000; % CHIRP bandwidth (Hz)

deltaMax = (delta+Pchirp)/2;        % little effective deltamax time in s
deltaMin = (delta-Pchirp)/2;        % little effective deltamin time in s

CHIRPtimeDelay = rampPct * Pchirp + preCHIRPdelay;

T = tD;                             % Sample time
Fs = 1/T;                           % Sampling frequency 
L = (nPts-omitPts)*(2^zf);          % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

echoVec = tE*(omitEchoes+1):tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/280.47;                       % um, 280.47 Hz/um (for PM25)

%% Import CHIRP data
[ap , spec] = readTecmag4d(strcat(datadir,datafile,'.tnt'));
CHIRPdat = reshape(spec, nPts, nEchoes);
CHIRPdat = CHIRPdat(1:end-omitPts,omitEchoes+1:end);

%% SNR calc 
n1 = CHIRPdat(1:noisePoints,:);
n2 = CHIRPdat(nPts-noisePoints-omitPts:end,:);
n = cat(1,n1,n2);
n = reshape(n,1,(2*noisePoints+1)*(nEchoes-omitEchoes));
s = reshape(CHIRPdat,1,(nPts-omitPts)*(nEchoes-omitEchoes));

figure
hold on
plot(abs(s))
plot(abs(n))


S = max(abs(s));
N = rms(n);

SNR = S/N
SNR_perRtScans = SNR/sqrt(2*ap.ns)

%% Apodization, zero filling, do FFT
pVec = 1:1:(nPts-omitPts);
filt = exp(-(pVec-(nPts-omitPts)/2).^2/((nPts-omitPts)/apofac)^2);
filt = repmat(filt',1,nEchoes-omitEchoes);

if apodize == 1
    CHIRPdat = CHIRPdat .* filt;
end

CHIRPdat = padarray(CHIRPdat, size(CHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

T2Dprofiles = fftshift(fft(CHIRPdat,NFFT)/L, 1); % Performs FFT algorithm

%% Plot CHIRP results
figure(1)
subplot(1,2,1)
plot(t*1e6,real(CHIRPdat(:,1)));
xlabel('time [us]')
subplot(1,2,2)
plot(z,2*abs(T2Dprofiles(:,1)),'LineWidth',1.5);
xlabel('real space [um]')
title('Plot of first T2-D FFT Profile and Echo')

figure(2)
surf(echoVec'/1000,z,abs(T2Dprofiles));
shading flat;
title('Surface plot of T2-D FFT Profiles')
xlabel('T2 [ms]')
ylabel('z [um]')
view([0 90])

%% No CHIRP load section
[ap,spec] = readTecmag4d(strcat(datadir,noCHIRPfile,'.tnt'));
data = reshape(spec,nPts,nEchoes);

% No CHIRP raw data and fft profiles
noCHIRPdat = reshape(data, nPts, nEchoes);
noCHIRPdat = noCHIRPdat(1:end-omitPts,omitEchoes+1:end);
if apodize == 1
    noCHIRPdat = noCHIRPdat .* filt;
end

noCHIRPdat = padarray(noCHIRPdat, size(noCHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's
CPprofiles = fftshift(fft(noCHIRPdat,NFFT)/L,1);

%% Plot first reference profile and coil profile

figure(3)
subplot(1,2,1)
plot(t*1e6,real(noCHIRPdat(:,1)));
xlabel('time [us]')
subplot(1,2,2)
plot(z,2*abs(CPprofiles(:,1)),'LineWidth',1.5);
xlabel('real space [um]')
title('Plot of first reference FFT Profile and Echo')

figure(4)
surf(echoVec'/1000,z,abs(CPprofiles));
shading flat;
title('Surface plot of reference FFT Profile')
xlabel('T2 [ms]')
ylabel('z [um]')
view([0 90])

figure(5)
hold on
plot(z,abs(CPprofiles(:,1))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','k')
plot(z,abs(T2Dprofiles(:,1))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','r')
hold off
xlabel('z [um]','fontsize',12)
title('T2-D and coil reference profiles')
set(gca,'Fontsize',12,'linewidth',2)

%% Coil Sensitivity Correction

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

%% Find Optimal data range with these figures

figure(7)
plot(abs(T2Dprofiles(:,1)))

t1_fig7=Pchirp*(BWchirp/2-f)/BWchirp;
deltaFig = 2*Pchirp*(BWchirp/2-f)/BWchirp + deltaMin; % expression for delta(effective), maybe

ylimits = [0 0.5];
deltaEff = 2*t1_fig7 ;
% deltaEff = delta - t1_fig7 - preCHIRPdelay;
deltaEff = fliplr(deltaEff);

figure(8)
subplot(4,1,1)
plot(abs(T2Dprofcorr(:,1)))
xlim([0 NFFT])
ylim(ylimits)
xlabel('index')

subplot(4,1,2)
plot(1e6*t1_fig7,abs(T2Dprofcorr(:,1)))
line(1e6*[0 0],[0 ylimits(2)])
line(1e6*[Pchirp Pchirp],[0 ylimits(2)])
line(1e6*[0+CHIRPtimeDelay 0+CHIRPtimeDelay],[0 ylimits(2)],'Color','r','LineStyle','--')
line(1e6*[Pchirp-CHIRPtimeDelay Pchirp-CHIRPtimeDelay],[0 ylimits(2)],'Color','r','LineStyle','--')
xlim(1e6*[min(t1_fig7), max(t1_fig7)]);
ylim(ylimits)
set(gca,'XDir','reverse')
xlabel('CHIRP time (us)')

subplot(4,1,3)
plot(1e6*deltaEff,abs(T2Dprofcorr(:,1)))
line(1e6*[0 0],[0 ylimits(2)])
line(1e6*[delta delta],[0 ylimits(2)])
line(1e6*[2*CHIRPtimeDelay 2*CHIRPtimeDelay],[0 ylimits(2)],'Color','r','LineStyle','--')
line(1e6*[delta-2*CHIRPtimeDelay delta-2*CHIRPtimeDelay ],[0 ylimits(2)],'Color','r','LineStyle','--')
xlim(1e6*[min(deltaEff), max(deltaEff)]);
ylim(ylimits)
% set(gca,'XDir','reverse')
xlabel('effective delta (us)')

subplot(4,1,4)
plot(z,abs(T2Dprofcorr(:,1)))
line([-sliceheight*1e3/2 -sliceheight*1e3/2],[0 ylimits(2)])
line([sliceheight*1e3/2 1e3*sliceheight/2],[0 ylimits(2)])
xlim([min(z), max(z)]);
ylim(ylimits)
% set(gca,'XDir','reverse')
xlabel('z (um)')
%% Data Range and Inversion

minind= 49;
maxind = 80;
% this is where I'm starting to put in some diffusion code. 

T2Ddat = abs(T2Dprofcorr(minind:maxind,:)); %crops data set according to above indices
deltaSteps = deltaEff(minind:maxind);

% deltaSteps = deltaFig(minind:maxind);

yD = log(T2Ddat./T2Ddat(1))';
xD = -gammaRad^2*G^2.*deltaSteps.^2.*(DELTA + (2/3)*deltaSteps);

T2Dsize = size(T2Ddat,1); % cuts down delta points to math those selected for the indices,assuming that the 

% fit to STE diffusion equation
p = polyfit(xD(1:T2Dsize),(yD(1,:)),1);

figure(9)
hold on
scatter(xD(1:T2Dsize),(yD(1,:))) 
plot(xD(1:T2Dsize),polyval(p,xD(1:T2Dsize)));

D = p(1)        % m^2 s^-1


%% surf of all T2-D Profiles

figure(10)
surf(echoVec/1000,deltaSteps*1e6,T2Ddat);
shading flat
colormap('jet');
colorbar 
xlabel('T2 [ms]'); 
ylabel('delta [us]');
title('D-T2 data')

%% Save data, display ILT Data params
close all

% T2Dexp = flipud(T2Ddat);
save(strcat(datadir,datafile, '.dat'), 'T2Ddat', '-ascii')

%UF Points [Min, Max; min(echoVec), max(echoVec), delta(eff)(min) [us], delta(eff)(max) [us], #echoes, #D points]
% sprintf('%f; %d %d %d; %.0f %.0f %.0f %.0f; %d %d',SNR, minind, maxind, firstinvertedind,  min(echoVec), max(echoVec), 1e6*min(t1), 1e6*max(t1), size(T1T2data,2), size(T1T2data,1))
dt = datestr(datetime('now','Format','dd MMMM yyyy HH:mm:ss'));
fileID = fopen(strcat(datadir,'DataNotesAuto.txt'),'a');
fprintf(fileID,'%s %s: %f; %d %d; %.0f %.0f %.2f %.2f; %d %d\n',dt,datafile, SNR, minind, maxind, min(echoVec), max(echoVec), 1e6*min(deltaSteps), 1e6*max(deltaSteps), size(T2Ddat,2), size(T2Ddat,1));
fclose(fileID);
