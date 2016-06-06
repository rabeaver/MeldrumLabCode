clear
clc
close all

%%
% CHIRP params
% ===================================
% ===== User-defined paramaters =====
% ===================================

datadir = 'Z:\VJL\';
datafile = 'CHIRP_DegassedWaterMolecularSieves_25mspw_sliceheight200um_tD8u_76pts_512scans_100nsWave_15April2016_result';
noCHIRPfile = 'noCHIRP_DegassedWaterMolecularSieves_25mspw_sliceheight200um_tD8u_76pts_512scans_100nsWave_15April2016_result';
filenameExt = '';

Pchirp = 0.025; % CHIRP Pulse Length (s)

sliceheight = 0.200; %mm
PreCPMGdelay = 100e-6; %s


nPts = 76; % # of acqu points
nEchoes = 32; % Echoes
tD = 8e-6; % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 700; %us
omitEchoes = 0; %the number of echoes to skip
noisePoints = 4; %number of points at beginning and end of each acqu period for noise
omitPts = 0; %blank spectrometer points to skip

zf = 1;                             % levels of zero filling
apodize = 0;                        %Gaussian apodization on (1) or off (0)?
apofac = 5;                         % Amount of Apodization

% ===================================
% === END User-defined paramaters ===
% ===================================

G = 6.59;                           % T m-1, B0 field gradient
gamma = 42.576;                     % MHz T-1
BWchirp = sliceheight*G*gamma*1000; % CHIRP bandwidth (Hz)

T = tD;                             % Sample time
Fs = 1/T;                           % Sampling frequency 
L = (nPts-omitPts)*(2^zf);      % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

echoVec = (omitEchoes+1)*tE:tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/280.47;                       % um, 280.47 Hz/um (for PM25)

%%
% Import CHIRP data
[ap , spec] = readTecmag4d(strcat(datadir,datafile,'.tnt'));

% CHIRPdat = spec(1,:);
% spec = spec2(nnn, :);
CHIRPdat = reshape(spec, nPts, nEchoes);
CHIRPdat = CHIRPdat(1:(end-omitPts),omitEchoes+1:end);

%% SNR calc (two sections)
n1 = CHIRPdat(1:noisePoints,:);
n2 = CHIRPdat(nPts-noisePoints-omitPts:end,:);
n = cat(1,n1,n2);
n = reshape(n,1,(2*noisePoints+1)*(nEchoes-omitEchoes));
s = reshape(CHIRPdat,1,(nPts-omitPts)*(nEchoes-omitEchoes));

figure
hold on
plot(abs(n))
plot(abs(s))

S = max(abs(s));
N = rms(n);

SNR = S/N
SNR_perRtScans = SNR/sqrt(2*ap.ns)
%%
pVec = 1:1:(nPts-omitEchoes);
filt = exp(-(pVec-(nPts-omitEchoes)/2).^2/((nPts-omitEchoes)/apofac)^2);
filt = repmat(filt',1,nEchoes);

if apodize == 1
    CHIRPdat = CHIRPdat .* filt;
end

CHIRPdat = padarray(CHIRPdat, size(CHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

T1T2profiles = fftshift(fft(CHIRPdat,NFFT)/L, 1); % Performs FFT algorithm

figure(1)
subplot(1,2,1)
hold on
plot(t*1e6,real(CHIRPdat(:,1)));
xlabel('time [us]')

subplot(1,2,2)
hold on
plot(z,2*abs(T1T2profiles(:,1)),'LineWidth',1.5);
xlabel('real space [um]')
title('Plot of first T1T2 FFT Profile and Echo')

hold off
figure(2)
hold on
surf(abs(T1T2profiles));
shading flat;
title('Surface plot of T1T2 FFT Profiles')
hold off

%% No CHIRP load section
close all

[~,spec] = readTecmag4d(strcat(datadir,noCHIRPfile,'.tnt'));
data = reshape(spec,nPts,nEchoes);

% No CHIRP raw data and fft profiles
% data = spec2(2,:);
noCHIRPdat = reshape(data, nPts, nEchoes);
noCHIRPdat = noCHIRPdat(1:(end-omitPts),omitEchoes+1:end);
if apodize == 1
    noCHIRPdat = noCHIRPdat .* filt;
end

noCHIRPdat = padarray(noCHIRPdat, size(noCHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

CPprofiles = fftshift(fft(noCHIRPdat,NFFT)/L,1);

%% Plot first T1-T2 profile and coil profile

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


%% Plot first T1T2 profile and coil profile

close all

figure(5)
hold on
plot(z,abs(CPprofiles(:,1))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','k')
plot(z,abs(T1T2profiles(:,1))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','r')
hold off

xlabel('{\it z} (um)','fontsize',12)
title('T1-T2 & coil profiles')
set(gca,'Fontsize',12,'linewidth',2)

%% Coil Sensitivity Correction

for k = 1:nEchoes-omitEchoes
    pcorr(:,k) = abs(CPprofiles(:,1));
end

T1T2profcorr = T1T2profiles./pcorr;

% figure(6)
% pcolor(abs(T1T2profcorr)); 
% colormap('jet');
% shading interp;
% colorbar('linewidth',2)
% caxis([0 1])
% title('Coil sensitivity corrected T1-T2 profiles')

%% Find Optimal data range with these figures
close all

figure(8)
plot(abs(T1T2profiles(:,1)))

t1_fig7=Pchirp*(BWchirp/2-f)/BWchirp;

figure(7)
subplot(2,1,1)
plot(abs(T1T2profcorr(:,1)))
xlim([0 NFFT])
ylim([0 1.05])
subplot(2,1,2)
plot(t1_fig7,abs(T1T2profcorr(:,1)))
line([0 0],[-2 2])
line([Pchirp Pchirp],[-2 2])
xlim([min(t1_fig7), max(t1_fig7)]);
ylim([0 1.05])
set(gca,'XDir','reverse')
xlabel('CHIRPtime (s)')

%% Data Range and Inversion

% manually select indices for data range and inversion (zero point)
minind= 54;
maxind = 87;

T1T2profiles2=zeros((maxind-minind+1),nEchoes-omitEchoes);

for ii = 1:nEchoes-omitEchoes;
    [~,firstinvertedind] = min(abs(T1T2profcorr(minind:maxind,ii)));
    firstinvertedind = firstinvertedind+minind-1;
    T1T2profiles2(1:firstinvertedind-minind+1,ii) = (abs(T1T2profcorr(minind:firstinvertedind,ii)));
    T1T2profiles2(firstinvertedind-minind+2:end,ii) = -(abs(T1T2profcorr(firstinvertedind+1:maxind,ii)));
end


% T1T2data=T1T2profiles2;
T1T2data=T1T2profiles2/max(max(abs(T1T2profiles2)));
t1=Pchirp*(BWchirp/2-f(minind:maxind))/BWchirp;

%plot first T1 column
figure
scatter(t1*1000,T1T2data(:,1),'linewidth',2)
xlabel('{\it t}_1 (ms)','fontsize',30)
title('T1-T2, first T1 column')
set(gca,'Fontsize',30,'linewidth',2)
xlim([0 1000*Pchirp])
ylim([-1.1 1.1])

%% surf of all T1-T2 Profiles

figure
surf(echoVec(:,1:end)*1000,t1*1000,T1T2data(:,1:end)); 
shading flat;
colormap('jet');
% shading interp;
colorbar 
ylabel('{\it t}_1 (ms)'); 
xlabel('{\it t}_2 (ms)');
title('T1-T2 data')


%% Save data, display ILT Data params
close all

T1T2data = T1T2data(:,1:end);
T1T2data2 = flipud(T1T2data);
data2d = T1T2data2;
save(strcat(datadir,datafile,filenameExt, '.dat'), 'T1T2data2', '-ascii');
% size(T1T2data);
% 1e6*abs(t1(1)-t1(end)); %#ok<NOPTS>
% 1e6*[min(t1), max(t1)]; %#ok<NOPTS>

%UF Points [Min, Max, Inv; min(echoVec), max(echoVec), InvTime(min) [us], InvTime(max) [us], #echoes, #T1 points]
% sprintf('%f; %d %d %d; %.0f %.0f %.0f %.0f; %d %d',SNR, minind, maxind, firstinvertedind,  min(echoVec), max(echoVec), 1e6*min(t1), 1e6*max(t1), size(T1T2data,2), size(T1T2data,1))

fileID = fopen(strcat(datadir,'DataNotesAuto.txt'),'a');
fprintf(fileID,'%s: %f; %d %d %d; %.0f %.0f %.0f %.0f; %d %d\n',datafile, SNR, minind, maxind, firstinvertedind,  min(echoVec), max(echoVec), 1e6*min(t1), 1e6*max(t1), size(T1T2data,2), size(T1T2data,1));
fclose(fileID);

