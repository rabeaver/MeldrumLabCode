
clear
clc
close all

%%
% CHIRP params
% ===================================
% ===== User-defined paramaters =====
% ===================================
%


spectrometer = 'Tecmag'; %'Tecmag' OR 'Kea'
datadir = 'C:\CommonData\TKM\Gouda\';
datafile = 'Gouda_CHIRP_2_13Apr2017_result'; %\1\data'; 
noCHIRPfile = 'Gouda_noCHIRP_1_13Apr2017_result'; %\1\data'; 


Pchirp = 1727e-6;                  % CHIRP Pulse Length (s)
pw     = 6e-6;                      % hard pulse length
sliceheight = 0.300;                % mm

rampPct = 0.0;                     % percent for the CHIRP power ramp to reach pMax

nPts = 36;                          % # of acqu points

omitPtsBack = 0;                    % the number of points at the end of each echo window that are zeros from the spectrometer
omitPtsFront = 0;                    % the number of points at the beginning of each echo window to zero
nEchoes = 512;                      % Echoes
omitEchoes = 0;                     % number of echoes to remove from data
echoChoice = 2;                     %the echo to use for display purposes

tD = 6e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 302-6;                           % us
preCHIRPdelay = 20e-6;             % s
noisePoints = 2;                    % number of points for measuring noise

nScans = 32768;                     % Number of scans in the experiment
cutRefPts = 0;                     %if necessary, can cut the data from the reference scan by half this value on each end of the acq window
                                    %use only if nPts for CHIRP on and CHIRP off expts don't match

zf = 2;                             % levels of zero filling

apodize = 0;                        % Gaussian apodization on (1) or off (0)?
apofac = 5;                         % Amount of Apodizatio


delta = 3.5e-3;                       % little delta time (s)
DELTA = 1e-3;                       % Big delta time in s


% ===================================
% === END User-defined paramaters ===
% ===================================

if strcmp(spectrometer,'Tecmag')==1
    G = 6.59;                           % T m-1, B0 field gradient
elseif strcmp(spectrometer,'Kea')==1    
    G = 23.87;
end

gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1
BWchirp = sliceheight*G*gamma*1000; % CHIRP bandwidth (Hz)

deltaMax = (delta+Pchirp)/2;        % little effective deltamax time in s
deltaMin = (delta-Pchirp)/2;        % little effective deltamin time in s

CHIRPtimeDelay = rampPct * Pchirp + preCHIRPdelay;

T = tD;                             % Sample time
Fs = 1/T;                           % Sampling frequency 
L = (nPts-omitPtsFront-omitPtsBack)*(2^zf);          % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

echoVec = tE*(omitEchoes+1):tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/(gamma*G);                    % um, 280.47 Hz/um (for PM25)

%% Import CHIRP data
if strcmp(spectrometer,'Tecmag')==1
    [ap , spec] = readTecmag4d(strcat(datadir,datafile,'.tnt'));
elseif strcmp(spectrometer,'Kea')==1
    [ap , spec] = readKea4d(strcat(datadir,datafile,'.2d'));
end

CHIRPdat = reshape(spec, nPts, nEchoes);
CHIRPdat = CHIRPdat(1+omitPtsFront:end-omitPtsBack,omitEchoes+1:end);

%% SNR calc 
n1 = CHIRPdat(1+omitPtsFront:noisePoints+omitPtsFront,:);
n2 = CHIRPdat(nPts-noisePoints-omitPtsBack-omitPtsFront:end,:);
n = cat(1,n1,n2);
n = reshape(n,1,(2*noisePoints+1)*(nEchoes-omitEchoes));
s = reshape(CHIRPdat,1,(nPts-omitPtsFront-omitPtsBack)*(nEchoes-omitEchoes));

figure
hold on
plot(abs(s))
plot(abs(n))


S = max(abs(s));
N = rms(n);

SNR = S/N
SNR_perRtScans = SNR/sqrt(nScans)

%% Apodization, zero filling, do FFT
pVec = 1:1:(nPts-omitPtsBack-omitPtsFront);
filt = exp(-(pVec-(nPts-omitPtsBack-omitPtsFront)/2).^2/((nPts-omitPtsBack-omitPtsFront)/apofac)^2);
filt = repmat(filt',1,nEchoes-omitEchoes);

if apodize == 1
    CHIRPdat = CHIRPdat .* filt;
end

CHIRPdat = padarray(CHIRPdat, size(CHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

T2Dprofiles = (fftshift(fft(CHIRPdat,NFFT)/L, 1)); % Performs FFT algorithm

%% Plot CHIRP results
figure(1)
subplot(1,2,1)
plot(t*1e6,real(CHIRPdat(:,echoChoice)));
xlabel('time [us]')
subplot(1,2,2)
% figure(5)
plot(z,2*abs(T2Dprofiles(:,echoChoice)),'LineWidth',1.5);
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
close all
if strcmp(spectrometer,'Tecmag')==1;
    [~ , spec] = readTecmag4d(strcat(datadir,noCHIRPfile,'.tnt'));
elseif strcmp(spectrometer,'Kea')==1;
    [~ , spec] = readKea4d(strcat(datadir,noCHIRPfile,'.2d'));
end

data = reshape(spec,nPts+cutRefPts,nEchoes);

% No CHIRP raw data and fft profiles
noCHIRPdat = reshape(data, nPts+cutRefPts, nEchoes);
noCHIRPdat = noCHIRPdat(1+omitPtsFront+cutRefPts/2:end-omitPtsBack-cutRefPts/2,omitEchoes+1:end);
if apodize == 1
    noCHIRPdat = noCHIRPdat .* filt;
end

noCHIRPdat = padarray(noCHIRPdat, size(noCHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's
CPprofiles = (fftshift(fft(noCHIRPdat,NFFT)/L,1));

%% Plot first reference profile and coil profile

figure(3)
subplot(1,2,1)
hold on
% plot(t*1e6,abs(noCHIRPdat(:,4)));
plot(t*1e6,real(noCHIRPdat(:,echoChoice)));
plot(t*1e6,imag(noCHIRPdat(:,echoChoice)));
xlabel('time [us]')
subplot(1,2,2)
plot(z,2*abs(CPprofiles(:,echoChoice)),'LineWidth',1.5);
xlabel('real space [um]')
title('Plot of first reference FFT Profile and Echo')

figure(4)
surf(echoVec'/1000,z,abs(CPprofiles));
shading flat;
title('Surface plot of reference FFT Profile')
xlabel('T2 [ms]')
ylabel('z [um]')
view([0 90])

hh = figure(5);
hold on
plot(z,abs(CPprofiles(:,echoChoice))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','k')
plot(z,abs(T2Dprofiles(:,echoChoice))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','r')
line([-sliceheight*1e3/2 -sliceheight*1e3/2],[0 2])
line([sliceheight*1e3/2 1e3*sliceheight/2],[0 2])
ylim([0 1.2])
hold off
xlabel('z [um]','fontsize',12)
title('T2-D and coil reference profiles')
set(gca,'Fontsize',12,'linewidth',2)
legend('ref','exp')
pubgraph(hh,14,2,'w','Arial')

%% Coil Sensitivity Correction

for k = 1:nEchoes-omitEchoes
    pcorr(:,k) = abs(CPprofiles(:,echoChoice));
end

T2Dprofcorr = T2Dprofiles./pcorr;

% figure(6)
% pcolor(abs(T1T2profcorr)); 
% colormap('jet');
% shading interp;
% colorbar('linewidth',2)
% caxis([0 1])
% title('Coil sensitivity corrected T1-T2 profiles')


% Calculate other axes; this version is for picking points automatically
% BigIndex = 1:NFFT;
% 
% FOV = 1/(gamma*1e6*G*tD);
% m_per_pt = FOV/NFFT;
% BigDELTA = DELTA + delta;
% 
% ptIndex = (1:NFFT);
% zBigIndex = FOV/2-BigIndex*m_per_pt;
% fBigIndex = zBigIndex * gamma*1e6 * G;
% 
% ptIndex = find(abs(fBigIndex)<=BWchirp/2);
% fIndex = zBigIndex(min(ptIndex):max(ptIndex)) * gamma*1e6 * G;
% deltaEffIndex = (1-(((BWchirp/2)-fIndex)/BWchirp))*2*Pchirp*1000;
% qIndex = 2*pi*gamma*1e6*G*deltaEffIndex/1000;
% vIndex = qIndex.^2.*(BigDELTA-deltaEffIndex./3000).*1e-9;
% T2Ddat = abs(T2Dprofcorr(minind:maxind,:)); %crops data set according to above indices

%% Find Optimal data range with these figures

close all 
figure(7)
plot(abs(T2Dprofiles(:,echoChoice)))

t1_fig7=Pchirp*(BWchirp/2-f)/BWchirp;
deltaFig = 2*Pchirp*(BWchirp/2-f)/BWchirp + deltaMin; % expression for delta(effective), maybe

wurstAmp = 1-(cos(pi*(t1_fig7)/Pchirp)).^40;


ylimits = [0 2];
deltaEff = 2*t1_fig7 ;
% deltaEff = delta - t1_fig7 - preCHIRPdelay;
deltaEff = fliplr(deltaEff);

figure(8)
subplot(4,1,1)
plot(abs(T2Dprofcorr(:,echoChoice)))
xlim([0 NFFT])
ylim(ylimits)
xlabel('index')
% 
subplot(4,1,2)
hold on
plot(1e6*t1_fig7,abs(T2Dprofcorr(:,echoChoice)))
plot(1e6*t1_fig7,wurstAmp,'--k');
line(1e6*[0 0],[0 ylimits(2)])
line(1e6*[Pchirp Pchirp],[0 ylimits(2)])
line(1e6*[0+CHIRPtimeDelay 0+CHIRPtimeDelay],[0 ylimits(2)],'Color','r','LineStyle','--')
line(1e6*[Pchirp-CHIRPtimeDelay Pchirp-CHIRPtimeDelay],[0 ylimits(2)],'Color','r','LineStyle','--')
xlim(1e6*[min(t1_fig7), max(t1_fig7)]);
ylim(ylimits)
set(gca,'XDir','reverse')
xlabel('CHIRP time (us)')

subplot(4,1,3)
plot(1e6*deltaEff,abs(T2Dprofcorr(:,echoChoice)))
line(1e6*[0 0],[0 ylimits(2)])
line(1e6*[delta delta],[0 ylimits(2)])
line(1e6*[2*CHIRPtimeDelay 2*CHIRPtimeDelay],[0 ylimits(2)],'Color','r','LineStyle','--')
line(1e6*[delta-2*CHIRPtimeDelay delta-2*CHIRPtimeDelay ],[0 ylimits(2)],'Color','r','LineStyle','--')
xlim(1e6*[min(deltaEff), max(deltaEff)]);
ylim(ylimits)
% set(gca,'XDir','reverse')
xlabel('effective delta (us)')

subplot(4,1,4)
plot(z,abs(T2Dprofcorr(:,echoChoice)))
line([-sliceheight*1e3/2 -sliceheight*1e3/2],[0 ylimits(2)])
line([sliceheight*1e3/2 1e3*sliceheight/2],[0 ylimits(2)])
xlim([min(z), max(z)]);
ylim(ylimits)
% set(gca,'XDir','reverse')
xlabel('z (um)')
%% Data Range and Inversion

minind = 66; %min(ptIndex);
maxind = 81; %max(ptIndex); 

% Calculate other axes
BigIndex = 1:NFFT;

FOV = 1/(gamma*1e6*G*tD);
m_per_pt = FOV/NFFT;
BigDELTA = DELTA + delta;

zBigIndex = FOV/2-BigIndex*m_per_pt;
fBigIndex = zBigIndex * gamma*1e6 * G;

ptIndex = minind:maxind; % find(abs(fBigIndex)<=BWchirp/2);
fIndex = zBigIndex(min(ptIndex):max(ptIndex)) * gamma*1e6 * G;
deltaEffIndex = (1-(((BWchirp/2)-fIndex)/BWchirp))*2*Pchirp*1000;
qIndex = 2*pi*gamma*1e6*G*deltaEffIndex/1000;
vIndex = qIndex.^2.*(BigDELTA-deltaEffIndex./3000).*1e-9;
T2Ddat = abs(T2Dprofcorr(minind:maxind,:)); %crops data set according to above indices



% deltaSteps = deltaEff(minind:maxind);

% deltaSteps = deltaFig(minind:maxind);

yD = log(T2Ddat./max(max((T2Ddat))))';
% xD = -gammaRad^2*G^2.*deltaSteps.^2.*(DELTA + (2/3)*deltaSteps);

T2Dsize = size(T2Ddat,1); % cuts down delta points to math those selected for the indices,assuming that the 

% fit to STE diffusion equation
p = polyfit(-vIndex, yD(1,:),1); %T2Dsize),(yD(1,:)),1);

figure(9)
hold on
scatter(-vIndex(1:T2Dsize),(yD(1,:)));
plot(-vIndex(1:T2Dsize),polyval(p,-vIndex(1:T2Dsize)));

D = p(1)*1e-9 ;        % *10-9 m^2 s^-1

%% surf of all T2-D Profiles

figure(10)
% surf(echoVec/1000,deltaSteps*1e6,T2Ddat);
surf(echoVec/1000,deltaEffIndex(1:end-5),T2Ddat(1:end-5,:));
shading flat
colormap('jet');
colorbar 
xlabel('T2 [ms]'); 
ylabel('delta [us]');
title('D-T2 data')

%% Save data, display ILT Data params
% close all

t2axis = echoVec; %s
% vaxis = gammaRad^2*G^2.*deltaSteps.^2.*((DELTA+delta) - (1/3)*deltaSteps); %s/m2
t2axis = t2axis';

% vaxis = [1.15, 1.39, 1.66, 1.94, 2.25, 2.58, 2.92, 3.29, 3.67, 4.08, 4.5, 4.93, 5.39, 5.86, 6.34]*1e7/1e9;

vIndex = rot90(vIndex,2)';


T2Ddat = ((T2Ddat));

save(strcat(datadir,datafile, '.dat'), 'T2Ddat', '-ascii')
save(strcat(datadir,datafile, '_T2axis.dat'), 't2axis', '-ascii')
save(strcat(datadir,datafile, '_vaxis.dat'), 'vIndex', '-ascii')

% %%
% 
% 
% Thmm = [0.001, 1]; %T2 (min and max)
% stepsh = 25; %horizontal steps
% Tvmm = [0.0001, 0.01]; %D min and max
% stepsv = 25;
% alpha = 1e7;
% 
% orient = 'b'; %both orientations
% kernel1 = 'exp(-h/T)';
% kernel2 = 'exp(-v*D)';
% 
% % Tvmm = Tvmm*1e9;
% 
% tic
%     [spectrum,tauh,tauv,chisq,compte]=upnnlsmooth3Dsvdfin(flipud(T2Ddat),echoVec*1e-6,rot90(vIndex,2),Thmm,stepsh,Tvmm,stepsv,alpha,-1,orient,kernel1,kernel2);
% toc
% 
% % the vertical (D) axis isn't calibrated correctly, and I'm not sure what
% % the reason is.
% 
% %     spectrum = flipdim(spectrum,1);
%     tauv = 1./tauv;
% %     tauv = tauv*1e-9;
%     tauv = flipdim(tauv,2);
% %     spectrum = spectrum';
% % figure
% % surf(echoVec*1e-6,rot90(vIndex,2),flipud(T2Ddat))
% % shading flat
% 
% %%
% % taulv = log10(tauv);
% stb = size(tauv);
% % taulh = log10(tauh);
% sta = size(tauh);
% tauv = tauv';
% 
% % spectrum = flipdim(spectrum,1);
% % tauv = 1./tauv;
% % tauv = flipdim(tauv,2);
%     
% surf(tauh,tauv,spectrum)    
% set(gca,'XScale','log','YScale','log')
% shading interp;
% %     %set(gcf,'Renderer','zbuffer');
% % axis([tauh(1),tauh(sta(2)),tauv(1),tauv(stb(2))]);
% view([0 90])
% xlabel('T_2 [s]')
% ylabel('D [m^2 s^{-1}]')
% 
% 
% 
% 
% %%
% 
% %UF Points [Min, Max; min(echoVec), max(echoVec), delta(eff)(min) [us], delta(eff)(max) [us], #echoes, #D points]
% % sprintf('%f; %d %d %d; %.0f %.0f %.0f %.0f; %d %d',SNR, minind, maxind, firstinvertedind,  min(echoVec), max(echoVec), 1e6*min(t1), 1e6*max(t1), size(T1T2data,2), size(T1T2data,1))
% dt = datestr(datetime('now','Format','dd MMMM yyyy HH:mm:ss'));
% fileID = fopen(strcat(datadir,'DataNotesAuto.txt'),'a');
% fprintf(fileID,'%s %s: %f; %d %d; %.0f %.0f %.2f %.2f; %d %d\n',dt,datafile, SNR, minind, maxind, min(echoVec), max(echoVec), 1e6*min(deltaSteps), 1e6*max(deltaSteps), size(T2Ddat,2), size(T2Ddat,1));
% fclose(fileID);


%%
omitEchoesT2D = 16;
omitDT2D = 5;
hlims = [0.001 1];
vlims = [0.001 1];
hsteps = 15;
vsteps = 15;
alpha = 1e7;

t2proc = t2axis(1+omitEchoesT2D:end)';
vproc = vIndex(1+omitDT2D:end)';
dataProc = T2Ddat(1+omitDT2D:end,1+omitEchoesT2D:end);

[spectrum,tauh,tauv,chisq,compte]=upnnlsmooth3Dsvdfin(dataProc,t2proc,vproc,hlims,hsteps,vlims,vsteps,alpha,-1,'b','exp(-h/T)','exp(-v*D)');

taulh = log10(tauh);
taulv = log10(tauv);

hh = figure(10);
surf(taulh,taulv,spectrum)
shading interp;
axis([taulh(1),taulh(end),taulv(1),taulv(end)]);
xmarks=(taulh(1):1:taulh(end));
ymarks=(taulv(1):1:taulv(end));

for aa = 0:range(xmarks)-2
    line([xmarks(aa+2) xmarks(aa+2)],[min(ymarks) max(ymarks)],[max(max(spectrum)) max(max(spectrum))],'LineStyle','--','LineWidth',2,'Color','w');
end

for aa = 0:range(ymarks)-2
    line([min(xmarks) max(xmarks)],[ymarks(aa+2) ymarks(aa+2)],[max(max(spectrum)) max(max(spectrum))],'LineStyle','--','LineWidth',2,'Color','w');
end


xlabel('$\log {\it{T}_2/\textrm{s}}$','FontSize',14,'interpreter','latex');
ylabel('$\log {\it{D}/  10^{-9} \textrm{m}^2 \textrm{s}^{-1}}$','FontSize',14,'interpreter','latex');
title('$\it{D}-\it{T}_2$ correlation','FontSize',16,'interpreter','latex');
pubgraph(hh,16,2,'w','mwa_cmr10')
