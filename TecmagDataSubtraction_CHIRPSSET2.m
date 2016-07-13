clear
clc
close all

%%
spectrometer = 'Tecmag'; %'Tecmag' OR 'Kea'
datadir = 'C:\Users\jnking01\Desktop\';
datafileWet = 'PureWater_tRepseries\Membrane_PureWater_CHIRP_7July2016_tRepSeries500';
datafileDry = 'Dry\Membrane_DryTreated_CHIRP_8July2016';
noCHIRPfileWet = 'PureWater_tRepseries\Membrane_PureWater_noCHIRP_7July2016_tRepseries';
noCHIRPfileDry = 'Dry\Membrane_DryTreated_noCHIRP_8July2016';



Pchirp = 196.8e-6;                  % CHIRP Pulse Length (s)
pw     = 6e-6;                      % hard pulse length
sliceheight = 0.200;                % mm
rampPct = 0.01;                     % percent for the CHIRP power ramp to reach pMax


nPts = 46;                          % # of acqu points
omitPtsBack = 0;                    % the number of points at the end of each echo window that are zeros from the spectrometer
omitPtsFront = 0;                    % the number of points at the beginning of each echo window to zero
nEchoes = 512;                      % Echoes
omitEchoes = 0;                     % numner of echoes to remove from data
tD = 2e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 180;                           % us
preCHIRPdelay = 0.2e-6;             % s
noisePoints = 5;                    % number of points for measuring noise

nScans = 16384;                      % Number of scans in the experiment
cutRefPts = 0;                     %if necessary, can cut the data from the reference scan by half this value on each end of the acq window
                                    %use only if nPts for CHIRP on and CHIRP off expts don't match

zf = 1;                             % levels of zero filling
apodize = 0;                        % Gaussian apodization on (1) or off (0)?
apofac = 5;                         % Amount of Apodizatio



delta = 0.40e-3;                       % little delta time (s)
DELTA = 0.0061e-3;                       % Big delta time in s

% ===================================
% === END User-defined paramaters ===
% ===================================

if strcmp(spectrometer,'Tecmag')==1;
    G = 6.59;                           % T m-1, B0 field gradient
elseif strcmp(spectrometer,'Kea')==1;    
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

%% Import CHIRP data from wet membrane
if strcmp(spectrometer,'Tecmag')==1;
    [ap , spec] = readTecmag4d(strcat(datadir,datafileWet,'.tnt'));
elseif strcmp(spectrometer,'Kea')==1;
    [ap , spec] = readKea4d(strcat(datadir,datafileWet,'.2d'));
end

CHIRPdatWet = reshape(spec, nPts, nEchoes);
CHIRPdatWet = CHIRPdatWet(1+omitPtsFront:end-omitPtsBack,omitEchoes+1:end);

%% Import CHIRP data from dry membrane
if strcmp(spectrometer,'Tecmag')==1;
    [ap , spec] = readTecmag4d(strcat(datadir,datafileDry,'.tnt'));
elseif strcmp(spectrometer,'Kea')==1;
    [ap , spec] = readKea4d(strcat(datadir,datafileDry,'.2d'));
end

CHIRPdatDry = reshape(spec, nPts, nEchoes);
CHIRPdatDry = CHIRPdatDry(1+omitPtsFront:end-omitPtsBack,omitEchoes+1:end);

%% Subtract Data and do FFT

CHIRPdatSub = CHIRPdatWet-CHIRPdatDry;

figure()
hold on
surf(abs(CHIRPdatSub))
shading flat

CHIRPdat = padarray(CHIRPdatSub, size(CHIRPdatSub(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

T2Dprofiles = flipud(fftshift(fft(CHIRPdat,NFFT)/L, 1)); % Performs FFT algorithm

%% Import noCHIRP data from wet membrane
close all
if strcmp(spectrometer,'Tecmag')==1;
    [~ , spec] = readTecmag4d(strcat(datadir,noCHIRPfileWet,'.tnt'));
elseif strcmp(spectrometer,'Kea')==1;
    [~ , spec] = readKea4d(strcat(datadir,noCHIRPfileWet,'.2d'));
end

data = reshape(spec,nPts+cutRefPts,nEchoes);

% No CHIRP raw data and fft profiles
noCHIRPdat = reshape(data, nPts+cutRefPts, nEchoes);
noCHIRPdatWet = noCHIRPdat(1+omitPtsFront+cutRefPts/2:end-omitPtsBack-cutRefPts/2,omitEchoes+1:end);

%% Import noCHIRP data from dry membrane
close all
if strcmp(spectrometer,'Tecmag')==1;
    [~ , spec] = readTecmag4d(strcat(datadir,noCHIRPfileDry,'.tnt'));
elseif strcmp(spectrometer,'Kea')==1;
    [~ , spec] = readKea4d(strcat(datadir,noCHIRPfileDry,'.2d'));
end

data = reshape(spec,nPts+cutRefPts,nEchoes);

% No CHIRP raw data and fft profiles
noCHIRPdat = reshape(data, nPts+cutRefPts, nEchoes);
noCHIRPdatDry = noCHIRPdat(1+omitPtsFront+cutRefPts/2:end-omitPtsBack-cutRefPts/2,omitEchoes+1:end);

%% Subtract data and perform FFT
noCHIRPdat = noCHIRPdatWet-noCHIRPdatDry;

noCHIRPdat = padarray(noCHIRPdat, size(noCHIRPdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's
CPprofiles = flipud(fftshift(fft(noCHIRPdat,NFFT)/L,1));

%% Plot first reference profile and coil profile

figure(3)
subplot(1,2,1)
plot(t*1e6,real(noCHIRPdat(:,4)));
xlabel('time [us]')
subplot(1,2,2)
plot(z,2*abs(CPprofiles(:,4)),'LineWidth',1.5);
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
plot(z,abs(CPprofiles(:,4))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','k')
plot(z,abs(T2Dprofiles(:,4))/max(abs(CPprofiles(:,1))),'linewidth',2,'color','r')
% ylim([0 0.1])
hold off
xlabel('z [um]','fontsize',12)
title('T2-D and coil reference profiles')
set(gca,'Fontsize',12,'linewidth',2)
legend('ref','exp')

%% Coil Sensitivity Correction

for k = 1:nEchoes-omitEchoes
    pcorr(:,k) = abs(CPprofiles(:,1));
end

T2Dprofcorr = T2Dprofiles./pcorr;

%%

% Need to only use points where abs(fIndex)<=BWchirp/2
BigIndex = 1:NFFT;

FOV = 1/(gamma*1e6*G*tD);
m_per_pt = FOV/NFFT;
BigDELTA = DELTA + delta;

ptIndex = (1:NFFT);
zBigIndex = FOV/2-BigIndex*m_per_pt;
fBigIndex = zBigIndex * gamma*1e6 * G;

ptIndex = find(abs(fBigIndex)<=BWchirp/2);
fIndex = zBigIndex(min(ptIndex):max(ptIndex)) * gamma*1e6 * G;
deltaEffIndex = (1-(((BWchirp/2)-fIndex)/BWchirp))*2*Pchirp*1000;
qIndex = 2*pi*gamma*1e6*G*deltaEffIndex/1000;
vIndex = qIndex.^2.*(BigDELTA-deltaEffIndex./3000).*1e-9;

%% Data Range and Inversion

minind = min(ptIndex);
maxind = max(ptIndex);
% this is where I'm starting to put in some diffusion code. 

T2Ddat = abs(T2Dprofcorr(minind:maxind,:)); %crops data set according to above indices
% deltaSteps = deltaEff(minind:maxind);

% deltaSteps = deltaFig(minind:maxind);

yD = log(T2Ddat./T2Ddat(1))';
% xD = -gammaRad^2*G^2.*deltaSteps.^2.*(DELTA + (2/3)*deltaSteps);

T2Dsize = size(T2Ddat,1); % cuts down delta points to math those selected for the indices,assuming that the 

% fit to STE diffusion equation
p = polyfit(-vIndex(1:T2Dsize),(yD(1,:)),1);

figure(9)
hold on
scatter(-vIndex(1:T2Dsize),(yD(1,:))) 
plot(-vIndex(1:T2Dsize),polyval(p,-vIndex(1:T2Dsize)));

D = p(1)        % *10-9 m^2 s^-1

%% surf of all T2-D Profiles

figure(10)
% surf(echoVec/1000,deltaSteps*1e6,T2Ddat);
surf(echoVec/1000,deltaEffIndex,T2Ddat);
shading flat
colormap('jet');
colorbar 
xlabel('T2 [ms]'); 
ylabel('delta [us]');
title('D-T2 data')

%% Save data, display ILT Data params
% close all

t2axis = echoVec*1e-6; %s
% vaxis = gammaRad^2*G^2.*deltaSteps.^2.*((DELTA+delta) - (1/3)*deltaSteps); %s/m2
t2axis = t2axis';

% vaxis = [1.15, 1.39, 1.66, 1.94, 2.25, 2.58, 2.92, 3.29, 3.67, 4.08, 4.5, 4.93, 5.39, 5.86, 6.34]*1e7/1e9;

vIndex = rot90(vIndex,2)';
vIndex = vIndex';
T2Ddat = flipud(T2Ddat);
% T2Dexp = flipud(T2Ddat);
save(strcat(datadir,datafileWet, '.dat'), 'T2Ddat', '-ascii')
save(strcat(datadir,datafileWet, '_T2axis.dat'), 't2axis', '-ascii')
save(strcat(datadir,datafileWet, '_vaxis.dat'), 'vIndex', '-ascii')