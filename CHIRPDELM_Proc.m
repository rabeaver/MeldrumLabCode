clear
clc
close all

%%
% CHIRP params
% ===================================
% ===== User-defined paramaters =====
% ===================================

Pchirp = 0.000622; % CHIRP Pulse Length (s)
sliceheight = 0.350; %mm

nPts = 76; % # of acqu points
nEchoes = 64; % Echoes
tD = 8e-6; % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 700; %us
omitEchoPts = 0; %the number of points that are zeros from the spectrometer
% nnn = 1; %expt number (for 2D CHIRP expts)

zf = 1;                             % levels of zero filling
apodize = 0;                        %Gaussian apodization on (1) or off (0)?
apofac = 5;                         % Amount of Apodization

vEcho = 1;                          %echo number for vieweing in plots

% ===================================
% === END User-defined paramaters ===
% ===================================

G = 6.59;                           % T m-1, B0 field gradient
gamma = 42.576;                     % MHz T-1
BWchirp = sliceheight*G*gamma*1000; % CHIRP bandwidth (Hz)

T = tD;                             % Sample time
Fs = 1/T;                           % Sampling frequency 
L = (nPts-omitEchoPts)*(2^zf);      % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

echoVec = tE:tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/280.47;                       % um, 280.47 Hz/um (for PM25)

%%
datadir = 'C:\Users\NMRLab\Desktop\CHIRP\DELM\';
datafile = 'CHIRP_DELM_Glycerol_test3_17Sept2015';

% Import CHIRP data
[~ , spec, spec2, ~] = readTecmag4d(strcat(datadir,datafile,'.tnt'));

% CHIRPdat = spec(1,:);
% spec = spec2(nnn, :);
CHIRPdat = reshape(spec, nPts, nEchoes);
CHIRPdat = CHIRPdat(1:end-omitEchoPts,:);

%%

pVec = 1:1:(nPts-omitEchoPts);
filt = exp(-(pVec-(nPts-omitEchoPts)/2).^2/((nPts-omitEchoPts)/apofac)^2);
filt = repmat(filt',1,nEchoes);

if apodize == 1
    CHIRPdat = CHIRPdat .* filt;
end

CHIRPdat = padarray(CHIRPdat, size(CHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

T1T2profiles = fftshift(fft(CHIRPdat,NFFT)/L, 1); % Performs FFT algorithm

figure(1)
subplot(1,2,1)
hold on
plot(t*1e6,real(CHIRPdat(:,vEcho)));
xlabel('time [us]')

subplot(1,2,2)
hold on
plot(z,2*abs(T1T2profiles(:,vEcho)),'LineWidth',1.5);
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

noCHIRPfile = 'noCHIRP_DELM_Glycerol_test3_17Sept2015';
[~,spec,spec2] = readTecmag4d(strcat(datadir,noCHIRPfile,'.tnt'));
data = reshape(spec,nPts,nEchoes);

% No CHIRP raw data and fft profiles
% data = spec2(2,:);
noCHIRPdat = reshape(data, nPts, nEchoes);
noCHIRPdat = noCHIRPdat(1:end-omitEchoPts,:);
if apodize == 1
    noCHIRPdat = noCHIRPdat .* filt;
end

noCHIRPdat = padarray(noCHIRPdat, size(noCHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's


CPprofiles = fftshift(fft(noCHIRPdat,NFFT)/L,1);

%% Plot first T1-T2 profile and coil profile

figure(3)
subplot(1,2,1)
hold on
plot(t*1e6,real(noCHIRPdat(:,vEcho)));
xlabel('time [us]')
subplot(1,2,2)
hold on
plot(z,2*abs(CPprofiles(:,vEcho)),'LineWidth',1.5);
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
plot(z,abs(CPprofiles(:,vEcho))/max(abs(CPprofiles(:,vEcho))),'linewidth',2,'color','k')
plot(z,abs(T1T2profiles(:,vEcho))/max(abs(CPprofiles(:,vEcho))),'linewidth',2,'color','r')
hold off

xlabel('{\it z} (um)','fontsize',12)
title('T1-T2 & coil profiles')
set(gca,'Fontsize',12,'linewidth',2)

%% Coil Sensitivity Correction

for k = 1:nEchoes
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
plot(abs(T1T2profiles(:,vEcho)))

t1_fig7=Pchirp*(BWchirp/2-f)/BWchirp;


figure(7)
subplot(2,1,1)
plot(abs(T1T2profcorr(:,vEcho)))
xlim([0 NFFT])
ylim([0 1.1])
subplot(2,1,2)
plot(t1_fig7,abs(T1T2profcorr(:,vEcho)))
line([0 0],[-2 2])
line([Pchirp Pchirp],[-2 2])
xlim([min(t1_fig7), max(t1_fig7)]);
ylim([0 1.1])
set(gca,'XDir','reverse')
xlabel('CHIRPtime (s)')



%% Data Range and Inversion

% manually select indices for data range and inversion (zero point)
minind= 100;
maxind = 150;
firstinvertedind = 110;

% automatically select indices
% minind=find(f>-BWchirp/2,1,'first');
% maxind=find(f<BWchirp/2,1,'last');
% [~,firstinvertedind] = min(abs(T1T2profiles(minind:maxind,3)));

T1T2profiles2=zeros((maxind-minind+1),nEchoes);
T1T2profiles2(1:firstinvertedind-minind+1,:) = (abs(T1T2profcorr(minind:firstinvertedind,:)));
T1T2profiles2(firstinvertedind-minind+2:end,:) = -(abs(T1T2profcorr(firstinvertedind+1:maxind,:)));

% T1T2data=T1T2profiles2;
T1T2data=T1T2profiles2/max(max(T1T2profiles2));
t1=Pchirp*(BWchirp/2-f(minind:maxind))/BWchirp;

%plot first T1 column
figure
scatter(t1*1000,T1T2data(:,vEcho),'linewidth',2)
xlabel('{\it t}_1 (ms)','fontsize',30)
title('T1-T2, first T1 column')
set(gca,'Fontsize',30,'linewidth',2)
% xlim([0 1000*Pchirp])
% ylim([-1.1 1.1])


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

%% T1 fit in cftool
echoNr = 1;
cftool(t1,T1T2data(:,echoNr));

%% Save data, display ILT Data params
close all

T1T2data = T1T2data(:,1:end);
T1T2data2 = flipud(T1T2data);
save(strcat(datadir,datafile, '.dat'), 'T1T2data2', '-ascii')
size(T1T2data)
1e6*abs(t1(1)-t1(end))
1e6*[min(t1), max(t1)]

%% T1Test
% For comparing your data to the data what you expect

close all

T1_1 = 0.0125; % T1 (s)
T1_2 = 0.0125;
w1 = 1; % Weights
w2 = 0;

t1new = linspace(max(t1), 0, length(t1)); % Simulated T1 Axis

% Make T1 Data
T1data1 = 1-2.*exp(-t1new./T1_1);
T1data2 = 1-2.*exp(-t1new./T1_2);

T1dat = w1.*T1data1 + w2.*T1data2;

figure()
hold on
plot(t1new, T1dat, '-r')
plot(t1, T1T2data(:,1), '*b')
hold off
