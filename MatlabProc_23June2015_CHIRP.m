clear
clc
close all

%%
% CHIRP params

Pchirp = 0.0020; % CHIRP Pulse Length (s)
BWchirp = 11223; % CHIRP bandwidth (Hz)

nPts = 69; % # of acqu points
nEchoes = 32; % Echoes
tD = 6e-6; % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 500; %us
omitEchoPts = 5; %the number of points that are zeros from the spectrometer
% nnn = 5; %expt number

zf = 2; % zero filling
T = tD*(2^zf);                     % Sample time
Fs = 1/T;                    % Sampling frequency
L = (nPts-omitEchoPts)*(2^zf);                     % Length of signal
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
apodize = 0; %Gaussian apodization on (1) or off (0)?


echoVec = tE:tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;                % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);          %Hz
z = f/280.47;           %um, 280.47 Hz/um (for PM25)

%%
datadir = '/Users/tyler/Desktop/';
datafile = 'CHIRP_GdH2O_5and50mM_20mspulse_20July2015_32768scans';

% Import CHIRP data
[~ , spec, spec2, ~] = readTecmag4d(strcat(datadir,datafile,'.tnt'));

% CHIRPdat = spec(1,:);
% spec = spec2(nnn, :);
CHIRPdat = reshape(spec, nPts, nEchoes);
CHIRPdat = CHIRPdat(1:end-omitEchoPts,:);

% figure(11)
% plot(abs(CHIRPdat(:,1)))

%%

pVec = 1:1:(nPts-omitEchoPts);
filt = exp(-(pVec-(nPts-omitEchoPts)/2).^2/((nPts-omitEchoPts)/5)^2);
% plot(filt)
filt = repmat(filt',1,nEchoes);
if apodize == 1
    CHIRPdat = CHIRPdat .* filt;
end
% surf(abs(CHIRPdat))

CHIRPdat = padarray(CHIRPdat, size(CHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

T1T2profiles = fftshift(fft(CHIRPdat,NFFT)/L, 1); % Performs FFT algorithm

figure(1)
subplot(1,2,1)
hold on
plot(t*1e6,real(CHIRPdat(:,2)));
% plot(t*1e6,imag(CHIRPdat(:,2)));
% plot(t*1e6,abs(CHIRPdat(:,2)));
xlabel('time [us]')

subplot(1,2,2)
hold on
plot(z,2*abs(T1T2profiles(:,2)),'LineWidth',1.5);
% plot(z,2*real(T1T2profiles(:,2)));
% plot(z,2*imag(T1T2profiles(:,2)));
xlabel('real space [um]')
title('Plot of Third T1T2 FFT Profile and Echo')

hold off
figure(2)
hold on
surf(abs(T1T2profiles));
shading flat;
title('Surface plot of T1T2 FFT Profiles')
hold off


%% No CHIRP load section
filenameNO = 'noCHIRP_GdH2O_5and50mM_20mspulse_20July2015_32768scans';
[~,spec,spec2] = readTecmag4d(strcat(datadir,filenameNO,'.tnt'));
data = reshape(spec,nPts,nEchoes);

% No CHIRP raw data and fft profiles
% data = spec2(1,:);
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
plot(t*1e6,real(noCHIRPdat(:,3)));
% plot(t*1e6,imag(noCHIRPdat(:,3)));
% plot(t*1e6,abs(noCHIRPdat(:,3)));
xlabel('time [us]')
subplot(1,2,2)
hold on
plot(z,2*abs(CPprofiles(:,3)),'LineWidth',1.5);
% plot(z,2*real(CPprofiles(:,3)));
% plot(z,2*imag(CPprofiles(:,3)));
xlabel('real space [um]')
title('Plot of Third T1T2 FFT Profile and Echo')
hold off

figure(4)
hold on
surf(abs(CPprofiles));
shading flat;
title('Surface plot of T1T2 FFT Profiles')
hold off


%% Plot first T1T2 profile and coil profile

figure(5)
hold on
plot(z,abs(CPprofiles(:,1))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','k')
plot(z,abs(T1T2profiles(:,1))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','r')
hold off

xlabel('{\it z} (um)','fontsize',12)
title('T1-T2 & coil profiles')
set(gca,'Fontsize',12,'linewidth',2)

%% Coil Sensitivity Correction

for k = 1:nEchoes
    pcorr(:,k) = abs(CPprofiles(:,1));
end

T1T2profcorr = T1T2profiles./pcorr;

figure(6)
pcolor(abs(T1T2profcorr)); 
colormap('jet');
shading interp;
set(gcf,'Renderer','painters');
colorbar('linewidth',2)
caxis([0 1])
title('Coil sensitivity corrected T1-T2 profiles')

%% Find Optimal data range with these figures
close all

figure(7)
plot(abs(T1T2profcorr(:,2)))

figure(8)
plot(abs(T1T2profiles(:,1)))

%% Data Range and Inversion

% manually select indices for data range and inversion (zero point)
minind= 114;
maxind = 140;
firstinvertedind = 132;

% automatically select indices
% minind=find(f>-BWchirp/2,1,'first');
% maxind=find(f<BWchirp/2,1,'last');
% [~,firstinvertedind] = min(abs(T1T2profiles(minind:maxind,3)));
% firstinvertedind = firstinvertedind + minind;


T1T2profiles2=zeros((maxind-minind+1),nEchoes);
T1T2profiles2(1:(firstinvertedind-minind),:)=abs(T1T2profcorr(minind:(firstinvertedind-1),:));
T1T2profiles2((firstinvertedind-minind+1):(maxind-minind+1),:)=-abs(T1T2profcorr(firstinvertedind:maxind,:))+repmat(abs(T1T2profcorr(firstinvertedind,:)), maxind-firstinvertedind+1, 1);
T1T2data=T1T2profiles2/max(max(T1T2profiles2));

t1=Pchirp*(BWchirp/2-f(minind:maxind))/BWchirp;

%plot first T1 column
figure
scatter(t1*1000,T1T2data(:,1),'linewidth',2)
xlabel('{\it t}_1 (ms)','fontsize',30)
title('T1-T2, first T1 column')
set(gca,'Fontsize',30,'linewidth',2)
% xlim([0 1000*Pchirp])
ylim([-1.1 1.1])

%% surf of all D-T2 Profiles

figure
surf(echoVec(:,1:end)*1000,t1*1000,T1T2data(:,1:end)); 
shading flat;
colormap('jet');
% shading interp;
%set(gcf,'Renderer','painters');
%daspect([1 1 1]);
%caxis([0 3e6])
colorbar %('linewidth',2)
ylabel('{\it t}_1 (ms)'); %,'fontsize',30)
xlabel('{\it t}_2 (ms)'); %,'fontsize',30)
title('T1-T2 data')
% set(gca,'Fontsize',30,'linewidth',2)
%set(gca,'XScale','log');
%set(gca,'YScale','log');
%% T1 fit
echoNr = 3;
cftool(t1,T1T2data(:,echoNr));
%%

T1T2data = T1T2data(:,1:end);
T1T2data2 = flipud(T1T2data);
save(strcat(datadir,datafile, '.dat'), 'T1T2data2', '-ascii')
size(T1T2data)
1e6*(t1(1)-t1(end))
1e6*[min(t1), max(t1)]