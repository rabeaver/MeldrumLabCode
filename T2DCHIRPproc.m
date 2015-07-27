clear
clc
close all

%%
% CHIRP params

Pchirp = 0.0025; % CHIRP Pulse Length (s)
BWchirp = 11223; % CHIRP bandwidth (Hz)

nPts = 69; % # of acqu points
nEchoes = 64; % Echoes
tD = 6e-6; % 2 * tD (Dwell time of 4e-06 should be input as 8e-06)
tE = 500; %us
omitEchoPts = 3; %the number of points that are zeros from the spectrometer
% nnn = 5; %expt number
eN = 3;

zf = 1; % zero filling
T = tD*(2^zf);                     % Sample time
Fs = 1/T;                    % Sampling frequency
L = (nPts-omitEchoPts)*(2^zf);                     % Length of signal
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
apodize = 0; %Gaussian apodization on (1) or off (0)?

gammaRad = 2.675222005e8; %s-1 T-1
gammaMHz = 42.57748; %MHz T-1
G = 6.59; % T m-1
DELTA = 1e-3; %s

echoVec = tE:tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;                % Time vector
f = linspace(Fs/2,-Fs/2,NFFT);          %Hz, flipped from [-Fs/2,+Fs/2] to [Fs/2, -Fs/2] to get alignment with our gradient, CHIRP frequency sweep direction. TKM, 7-17-2015
z = f/280.47;           %um, 280.47 Hz/um (for PM25)

%%
datadir = '/Users/jaredking/Documents/Chemistry/Research/CHIRP/';
datafile = 'T2D_STE_CHIRP_5dB_2048_Glycerol_2_17July2015';

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
filt = exp(-(pVec-(nPts-omitEchoPts)/2).^2/((nPts-omitEchoPts)/10)^2);
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
plot(t*1e6,real(CHIRPdat(:,eN)));
% plot(t*1e6,imag(CHIRPdat(:,2)));
% plot(t*1e6,abs(CHIRPdat(:,2)));
xlabel('time [us]')

subplot(1,2,2)
hold on
plot(z,2*abs(T1T2profiles(:,eN)),'LineWidth',1.5);
% plot(z,2*real(T1T2profiles(:,2)));
% plot(z,2*imag(T1T2profiles(:,2)));
xlabel('real space [um]')
title('Plot of Third T2D FFT Profile and Echo')

hold off
figure(2)
hold on
surf(abs(T1T2profiles));
shading flat;
title('Surface plot of T2D FFT Profiles')
hold off


%% No CHIRP load section
filenameNO = 'T2D_STE_noCHIRP_20dB_2048_Glycerol_2_17July2015';
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
plot(t*1e6,real(noCHIRPdat(:,eN)));
% plot(t*1e6,imag(noCHIRPdat(:,3)));
% plot(t*1e6,abs(noCHIRPdat(:,3)));
xlabel('time [us]')
subplot(1,2,2)
hold on
plot(z,2*abs(CPprofiles(:,eN)),'LineWidth',1.5);
% plot(z,2*real(CPprofiles(:,3)));
% plot(z,2*imag(CPprofiles(:,3)));
xlabel('real space [um]')
title('Plot of Third T2D FFT Profile and Echo')
hold off

figure(4)
hold on
surf(abs(CPprofiles));
shading flat;
title('Surface plot of T2D FFT Profiles')
hold off


%% Plot first T1T2 profile and coil profile

figure(5)
hold on
plot(z,abs(CPprofiles(:,eN))/max(abs(CPprofiles(:,eN))),'linewidth',2,'color','k')
plot(z,abs(T1T2profiles(:,eN))/max(abs(CPprofiles(:,eN))),'linewidth',2,'color','r')
hold off

xlabel('{\it z} (um)','fontsize',12)
title('T2-D & coil profiles')
set(gca,'Fontsize',12,'linewidth',2)

%% Coil Sensitivity Correction

for k = 1:nEchoes
    pcorr(:,k) = abs(CPprofiles(:,eN));
end

T1T2profcorr = T1T2profiles./pcorr;

figure(6)
pcolor(abs(T1T2profcorr)); 
colormap('jet');
shading interp;
set(gcf,'Renderer','painters');
colorbar('linewidth',2)
caxis([0 1])
title('Coil sensitivity corrected T2-D profiles')

%% Find Optimal data range with these figures
close all

figure(7)
plot(abs(T1T2profcorr(:,eN)))

figure(8)
plot(abs(T1T2profiles(:,eN)))

%% Data Range and Inversion

% manually select indices for data range and inversion (zero point)
minind= 114;
maxind = 124;
% firstinvertedind = 132;

% automatically select indices
% minind=find(f>BWchirp/2,1,'last');
% maxind=find(f<-BWchirp/2,1,'first');
% [~,firstinvertedind] = min(abs(T1T2profiles(minind:maxind,3)));
% firstinvertedind = firstinvertedind + minind;


% T1T2profiles2=zeros((maxind-minind+1),nEchoes);
% T1T2profiles2(1:(firstinvertedind-minind),:)=abs(T1T2profcorr(minind:(firstinvertedind-1),:));
% T1T2profiles2((firstinvertedind-minind+1):(maxind-minind+1),:)=-abs(T1T2profcorr(firstinvertedind:maxind,:))+repmat(abs(T1T2profcorr(firstinvertedind,:)), maxind-firstinvertedind+1, 1);
% T1T2data=T1T2profiles2/max(max(T1T2profiles2));

Ddat = abs(T1T2profcorr(minind:maxind,:));
Dvec=Pchirp*(BWchirp/2-f(minind:maxind))/BWchirp;

%plot first T1 column
figure
plot(Dvec*1000,Ddat(:,eN)','linewidth',2)
xlabel('{\it t}_1 (ms)','fontsize',30)
title('T2D')
ylabel(strcat('signal, echo ',num2str(eN),' [arb]')); %,'fontsize',30)
xlabel('delta [ms]');
set(gca,'Fontsize',30,'linewidth',2)
% xlim([0 1000*Pchirp])
% ylim([-1.1 1.1])

%% surf of all D-T2 Profiles

figure
surf(echoVec(:,eN:end)*1000,Dvec*1000,Ddat(:,eN:end)); 
shading flat;
colormap('jet');
% shading interp;
%set(gcf,'Renderer','painters');
%daspect([1 1 1]);
%caxis([0 3e6])
colorbar %('linewidth',2)
ylabel('delta (ms)'); %,'fontsize',30)
xlabel('{\it t}_2 (ms)'); %,'fontsize',30)
title('T2-D data')
% set(gca,'Fontsize',30,'linewidth',2)
%set(gca,'XScale','log');
%set(gca,'YScale','log');

figure
scatter(Dvec*1000,Ddat(:,eN))
ylabel(strcat('signal, echo ',num2str(eN),' [arb]')); %,'fontsize',30)
xlabel('delta [ms]');

%% T1 fit
% cftool(Dvec*1000,Ddat(:,eN));

%% Test for D from STE equation
yvals = log(Ddat(:,eN)./Ddat(1,eN));
xvals = -gammaRad^2*G^2*Dvec.^2.*(DELTA + (2/3)*Dvec);
% cftool(xvals,yvals)

% Given the followign diffusion coefficient, determine the correct value
% for the x axis
D = 3.5e-12; %m2 s-1

test = -yvals/(-gammaRad^2*G^2*D);

p3 = ones(length(test),1)*(2/3);
p2 = ones(length(test),1)*(1e-3);
p1 = ones(length(test),1)*(0);
p0 = test;

p = [p3,p2,p1,p0];

for i = 1:length(test)
    r(:,i) = roots(p(i,:));
end

deltaCalc = r(3,:)';
deltaVec = linspace(max(0,deltaCalc(2)),deltaCalc(end),length(test));

% retest the diffusion coefficient using this newly calculated vector
xvals2 = -gammaRad^2*G^2*deltaVec.^2.*(DELTA + (2/3)*deltaVec);
cftool(xvals2,yvals)
%%

Ddat = Ddat(:,eN:end);
T1T2data2 = (Ddat);
save(strcat(datafile, '.dat'), 'T1T2data2', '-ascii')
size(Ddat)
1e6*(Dvec(1)-Dvec(end))
1e6*[min(Dvec), max(Dvec)]