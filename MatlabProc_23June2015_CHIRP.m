clear
clc
close all

%%
% CHIRP params
%===================================
%===== User-defined paramaters =====
%===================================

Pchirp = 0.10; % CHIRP Pulse Length (s)
sliceheight = 0.300; %mm

nPts = 76; % # of acqu points
nEchoes = 128; % Echoes
tD = 8e-6; % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 700; %us
omitEchoPts = 2; %the number of points that are zeros from the spectrometer
% nnn = 1; %expt number (for 2D CHIRP expts)

zf = 1; % levels of zero filling
apodize = 0; %Gaussian apodization on (1) or off (0)?

%===================================
%=== END User-defined paramaters ===
%===================================

G = 6.59; %T m-1, B0 field gradient
gamma = 42.576; %MHz T-1
BWchirp = sliceheight*G*gamma*1000; % CHIRP bandwidth (Hz)

T = tD; %*(2^zf);                     % Sample time
Fs = 1/T;                    % Sampling frequency
L = (nPts-omitEchoPts)*(2^zf);                     % Length of signal
NFFT = 2^nextpow2(L); % Next power of 2 from length of y

echoVec = tE:tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;                % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);          %Hz
z = f/280.47;           %um, 280.47 Hz/um (for PM25)

% range(f)

%%
datadir = '/Users/jaredking/Documents/Chemistry/Research/CHIRP/';
datafile = 'CHIRP_Glycerol_40mspw_sliceheight20um_5db_25July2015';

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
filt = exp(-(pVec-(nPts-omitEchoPts)/2).^2/((nPts-omitEchoPts)/3)^2);
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
filenameNO = 'noCHIRP_Glycerol_40mspw_sliceheight20um_5db_25July2015';
[~,spec,spec2] = readTecmag4d(strcat(datadir,filenameNO,'.tnt'));
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
colorbar('linewidth',2)
caxis([0 1])
title('Coil sensitivity corrected T1-T2 profiles')

%% Find Optimal data range with these figures
close all

figure(8)
plot(abs(T1T2profiles(:,1)))

t1_fig7=Pchirp*(BWchirp/2-f)/BWchirp;


figure(7)
subplot(2,1,1)
plot(abs(T1T2profcorr(:,2)))
<<<<<<< HEAD
xlim([0 length(T1T2profiles(:,1))])
=======
xlim([0 NFFT])
ylim([0 1.1])
>>>>>>> e8a759cc27023e0b4278147a493c0f58d169c66f
subplot(2,1,2)
plot(t1_fig7,abs(T1T2profcorr(:,2)))
line([0 0],[-2 2])
line([Pchirp Pchirp],[-2 2])
xlim([min(t1_fig7), max(t1_fig7)]);
ylim([0 1.1])
set(gca,'XDir','reverse')
xlabel('CHIRPtime (s)')



%% Data Range and Inversion

% manually select indices for data range and inversion (zero point)
minind= 338;
maxind = 1652;
firstinvertedind = 1555;

% automatically select indices
% minind=find(f>-BWchirp/2,1,'first');
% maxind=find(f<BWchirp/2,1,'last');
% [~,firstinvertedind] = min(abs(T1T2profiles(minind:maxind,3)));

% firstinvertedind = firstinvertedind + minind;
% % firstinvertedind = NFFT/2;

T1T2profiles2=zeros((maxind-minind+1),nEchoes);
% T1T2profiles2((maxind-firstinvertedind+1):(maxind-minind+1),:) = (-abs(T1T2profcorr(minind:firstinvertedind,:))+repmat(abs(T1T2profcorr(firstinvertedind,:)), firstinvertedind-minind+1, 1));
% T1T2profiles2(1:(maxind-firstinvertedind),:) = (abs(T1T2profcorr((firstinvertedind+1):maxind,:)));
T1T2profiles2(1:firstinvertedind-minind+1,:) = (abs(T1T2profcorr(minind:firstinvertedind,:)));
T1T2profiles2(firstinvertedind-minind+2:end,:) = -(abs(T1T2profcorr(firstinvertedind+1:maxind,:)));

close all
T1guess = 0.052;
T1T2profilesTest = zeros(NFFT,1);
T1T2profilesTest(1:NFFT/2) = (abs(T1T2profcorr(1:NFFT/2)));
T1T2profilesTest(NFFT/2+1:NFFT) = -(abs(T1T2profcorr(NFFT/2+1:NFFT)));

% T1test = (1-2*exp((t1_fig7-Pchirp)/T1guess));
% T1test2 = (1-2*exp(-t1_fig7/T1guess));
% figure
% subplot(2,1,1)
% plot(T1T2profilesTest(:,1))
% xlim([0 NFFT])
% ylim([-2 2]);
% subplot(2,1,2)
% hold on
% plot(t1_fig7,T1T2profilesTest(:,1),'-k')
% % plot(t1_fig7,T1test,'-r')
% % plot(t1_fig7,T1test2,'-b')
% line([0 0],[-2 2])
% line([Pchirp Pchirp],[-2 2])
% % xlim([0, Pchirp]);
% % xlim([min(t1_fig7), 0.02]);
% ylim([-2 2]);
% xlabel('CHIRPtime (s)')
% set(gca,'XDir','reverse')

%


% T1T2profiles2(1:(firstinvertedind-minind),:) = abs(T1T2profcorr(minind:(firstinvertedind-1),:));
% T1T2profiles2((firstinvertedind-minind+1):(maxind-minind+1),:) = -abs(T1T2profcorr(firstinvertedind:maxind,:))+repmat(abs(T1T2profcorr(firstinvertedind,:)), maxind-firstinvertedind+1, 1);
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

%% Only using first half of data
minind=find(f>-BWchirp/2,1,'first');
[~,firstinvertedind] = min(abs(T1T2profiles(minind:maxind,3)));
firstinvertedind = firstinvertedind + minind;
T1T2data2 = (abs(T1T2profcorr(1:firstinvertedind,:)));
t2=Pchirp*(BWchirp/2-f(1:firstinvertedind))/BWchirp;
%plot first T1 column
figure
scatter(t2*1000,T1T2data2(:,1),'linewidth',2)
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


%% Normalize T2?
close all

T2_1 = 0.003; %T2 (s)
T2_2 = 0.0062;
echoVec2 = echoVec./1e6;

w1 = 0.9; % Weights
w2 = 0.3;

T2dat1 = w1.*exp(-echoVec2./T2_1);
T2dat2 = w2.*exp(-echoVec2./T2_2);

T2dat = T2dat1 + T2dat2;

figure()
plot( T2dat);

for i = 1:nEchoes
    T1T2datadiv(:,i) = T1T2data(:,i)/T2dat(i);
end

figure()
surf(T1T2datadiv); shading flat

%% T1Test

T1_1 = 0.0215; % T1 (s)
T1_2 = 0.0125;

w1 = 1; % Weights
w2 = 0;

t1eh = linspace(max(t1), 0, length(t1));

T1data1 = 1-2.*exp(-t1eh./T1_1);
T1data2 = 1-2.*exp(-t1eh./T1_2);

T1dat = w1.*T1data1 + w2.*T1data2;

% t1new=2*Pchirp*(BWchirp/2-f(minind:maxind))/BWchirp;
t1new = linspace(max(t1), 0, length(t1));


figure()
hold on
plot(t1eh, T1dat+1-.5857, '-r')
plot(t1, T1T2data(:,1), '*b')
hold off


%% T1 fit
echoNr = 3;
cftool(t1,T1T2data(:,echoNr));
%%

T1T2data = T1T2data(:,1:end);
T1T2data2 = flipud(T1T2data);
save(strcat(datadir,datafile, '.dat'), 'T1T2data2', '-ascii')
size(T1T2data)
1e6*abs(t1(1)-t1(end))
1e6*[min(t1), max(t1)]