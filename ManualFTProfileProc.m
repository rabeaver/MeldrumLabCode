clear
close all
clc
%%
% This script takes a series of CPMG files (from the CPMG sequence) at
% different lift positions, as well as a reference file (for coil
% sensititvity correction). It Fourier transforms each, does a sensitivity
% correction, then concatenatates the whole set into a 2D profile, but
% using FT rather than discrete measurement points.

% =====================================
% === BEGIN User-defined paramaters ===
% =====================================

spectrometer = 'Kea'; %'Tecmag'

%location of data sets
datadir = 'Z:\Data\TKM\GCI_Acrylate\17May2016\Acrylate_ManualProfile_FT\';
filenum = '1';
datafile = '\data';

%location of reference file
reffile = 'Z:\Data\TKM\GCI_Acrylate\17May2016\Rubber_FTRef\1\data';

%acq parameters

pw = 2e-6;                          % hard pulse length
nPts = 16;                          % # of acqu points
omitPtsBack = 0;                    % the number of points at the end of each echo window that are zeros from the spectrometer
omitPtsFront = 0;                   % the number of points at the beginning of each echo window to zero
nEchoes = 256;                      % Echoes
omitEchoes = 0;                     % numner of echoes to remove from data
tD = 2e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 62;                            % us
noisePoints = 4;                    % number of points for measuring noise
cutRefPts = 0;                      % if necessary, can cut the data from the reference scan by half this value on each end of the acq window
                                    % use only if nPts for CHIRP on and CHIRP off expts don't match
stepSize = 100;                     % um, lift motion increments
n2DPoints = 9;                      % number of manual lift positions

zf = 1;                             % levels of zero filling
apodize = 0;                        % Gaussian apodization on (1) or off (0)?
apofac = 5;                         % Amount of Apodization

fitopts = statset('MaxIter',5000,'TolX',1e-14,'UseParallel',true,'Display','off');

% ===================================
% === END User-defined paramaters ===
% ===================================

%set appropriate field gradients
if strcmp(spectrometer,'Tecmag')==1;
    G = 6.59;                           % T m-1, B0 field gradient
elseif strcmp(spectrometer,'Kea')==1;    
    G = 23.87;
end

gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1

%calculate FT frequency and position axes
T = tD;                             % Sample time
Fs = 1/T;                           % Sampling frequency 
L = (nPts-omitPtsFront-omitPtsBack)*(2^zf);          % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

echoVec = tE*(omitEchoes+1):tE:(nEchoes*tE);
t = (-(L-1)/2:L/2)*T;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/(gamma*G);                    % um, 280.47 Hz/um (for PM25)


%based on lift step size, calculate the number of points from each
%individual profile to concatenate into the superprofile
ftspacing = diff(z,1);
ftspacing = ftspacing(1,1);
nCenterPts = ceil(stepSize/ftspacing);    %number of points from each FT profile to use in the concatented superprofile.
ptRange = round([NFFT/2 - (nCenterPts-1)/2, NFFT/2 + (nCenterPts-1)/2]);
%% Import ref data

if strcmp(spectrometer,'Tecmag')==1;
        [refap , refspec] = readTecmag4d(strcat(reffile,'.tnt'));
    elseif strcmp(spectrometer,'Kea')==1;
        [refap , refspec] = readKea4d(strcat(reffile,'.2d'));
end

refdat = reshape(refspec, nPts, nEchoes);
refdat = refdat(1+omitPtsFront:end-omitPtsBack,omitEchoes+1:end);
    
% Apodization, zero filling, do FFT
pVec = 1:1:(nPts-omitPtsBack-omitPtsFront);
filt = exp(-(pVec-(nPts-omitPtsBack-omitPtsFront)/2).^2/((nPts-omitPtsBack-omitPtsFront)/apofac)^2);
filt = repmat(filt',1,nEchoes-omitEchoes);
    
if apodize == 1
    refdat = refdat .* filt;
end
    
refdat = padarray(refdat, size(refdat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's for zero filling
    
refprofile = flipud(fftshift(fft(refdat,NFFT)/L, 1)); % Performs FFT algorithm
refz = f/(gamma*G);

figure(1)
surf(echoVec,refz,abs(refprofile)); shading interp;

%make a matrix for sensitivity correction. This is simply the first echo of
%the reference experiment repeated to fill the appropriate size of the
%measured data
refcorr = repmat(abs(refprofile(:,1)),[1 nEchoes]);

%% Import measured data

%loop over each measured profile point
for ii = 1:n2DPoints
    if strcmp(spectrometer,'Tecmag')==1;
        [ap , spec] = readTecmag4d(strcat(datadir,num2str(ii),datafile,'.tnt'));
    elseif strcmp(spectrometer,'Kea')==1;
        [ap , spec] = readKea4d(strcat(datadir,num2str(ii),datafile,'.2d'));
    end
    
    dat = reshape(spec, nPts, nEchoes);
    dat = dat(1+omitPtsFront:end-omitPtsBack,omitEchoes+1:end);
    
    % Apodization, zero filling, do FFT
    pVec = 1:1:(nPts-omitPtsBack-omitPtsFront);
    filt = exp(-(pVec-(nPts-omitPtsBack-omitPtsFront)/2).^2/((nPts-omitPtsBack-omitPtsFront)/apofac)^2);
    filt = repmat(filt',1,nEchoes-omitEchoes);
    
    if apodize == 1
        dat = dat .* filt;
    end
    
    dat = padarray(dat, size(dat(:,1),1)/2*((2^zf)-1),0); % Pad with 0's
    
    profile(:,:,ii) = flipud(fftshift(fft(dat,NFFT)/L, 1))./refcorr; % Performs FFT algorithm and does sensitivity correction
    
    %I don't understand why the flipud is necessary, but it gives correct
    %results. I think it's just how the data is stored/manipulated.
    
    z(ii,:) = f/(gamma*G)-(ii-1)*(stepSize);                    % um, including the offset appropriate for the lift motion
    
end
    z = z';


%%  Plot individual lift positions in a series of subplots.
figure(2)
for jj = 1:n2DPoints
    subplot(ceil(sqrt(n2DPoints)),ceil(sqrt(n2DPoints)),jj)
    pcolor(echoVec*1e-3,z(ptRange(1):ptRange(2),jj),abs(profile(ptRange(1):ptRange(2),:,jj))); shading flat
    caxis([0 max(max(max(abs(profile(ptRange(1):ptRange(2),:,:)))))])
    ylim([round(min(z(:,end))) round(max(z(:,1)))]) %may need to figure out how to adjust the ylimits based on the step size and number of points
    xlim([0 nEchoes*tE*1e-3])
end
    xlabel('time [ms]')
    ylabel('position [um]')

%% Take portions of each profile for concatenation
%portions of the position axis
limZ = z(ptRange(1):ptRange(2),:);
limZ = flipud(limZ);
limZ = reshape(limZ,nCenterPts*n2DPoints,1);

%portions of each data set
limProfile = permute(profile,[1,3,2]);
limProfile = limProfile(ptRange(1):ptRange(2),:,:);
limProfile = flipud(limProfile);
limProfile = reshape(limProfile,nCenterPts*n2DPoints,nEchoes);

%plot all
figure(3)
surf(echoVec*1e-3,limZ,abs(limProfile)); shading flat
xlabel('time [ms]');
ylabel('position [um]');
view([0 90]);


%% Fit superprofile to exponential function(s)
guess = [1; 1];


% Monoexponential fit
for i=1:length(limZ)
    
   try
       [fit.beta(:,i),fit.residual(:,i),fit.J(:,:,i)] = nlinfit(echoVec*1e-3,abs(limProfile(i,:)),@t2monofit_simple,guess,fitopts);
%                final_data(k).beta(3,i) = 0;
%                final_data(k).beta(4,i) = 0;
   catch pm
       fit.beta(:,i) = [0;0];%;0;0];
       fit.J(:,:,i) = zeros(size(limProfile,1),2);
       fit.residual(:,i) = zeros(size(limProfile,1),1);
   end

   %guess(:,i) = final_data.beta(:,i);
   
   [fit.ypred(:,i),fit.delta(:,i)] = nlpredci(@t2monofit_simple,echoVec*1e-3,fit.beta(:,i),fit.residual(:,i),'Jacobian',fit.J(:,:,i));
   tt = nlparci(fit.beta(:,i),fit.residual(:,i),'jacobian',fit.J(:,:,i), 'alpha', 0.1);
   fit.ci(:,i) = tt(:,2)-fit.beta(:,i);
   [fit.yvals(:,i),fit.yvalsdelta(:,i)] = nlpredci(@t2monofit_simple,echoVec*1e-3,fit.beta(:,i),fit.residual(:,i),'Jacobian',fit.J(:,:,i));
end

%% Plot together

figure(4)
subplot(2,1,1)
hold on
plot(limZ,fit.beta(1,:),'-ok')
plot(limZ,fit.beta(1,:)+fit.ci(1,:),'--k')
plot(limZ,fit.beta(1,:)-fit.ci(1,:),'--k')
set(gca,'XDir','reverse')
set(gca,'YScale','log')
xlim([round(min(limZ)) round(max(limZ))])
ylim([0 1])
ylabel('signal amplitude [arb]')
subplot(2,1,2)
hold on
plot(limZ,fit.beta(2,:),'-or')
plot(limZ,fit.beta(2,:)+fit.ci(2,:),'--r')
plot(limZ,fit.beta(2,:)-fit.ci(2,:),'--r')
ylim([0 50])
xlim([round(min(limZ)) round(max(limZ))])
set(gca,'XDir','reverse')
set(gca,'YScale','log')
xlabel('position [um], top \leftarrow \rightarrow bottom')
ylabel('{\itT}_2 [ms]')


%%
figure(5)
subplot(2,1,1)
hold on
plot(limZ,fit.beta(1,:),'-ok')
plot(limZ,fit.beta(1,:)+fit.ci(1,:),'--k')
plot(limZ,fit.beta(1,:)-fit.ci(1,:),'--k')
set(gca,'XDir','reverse')
set(gca,'YScale','log')
xlim([-700 -400 ])
ylim([0 0.1])
ylabel('signal amplitude [arb]')
subplot(2,1,2)
hold on
plot(limZ,fit.beta(2,:),'-or')
plot(limZ,fit.beta(2,:)+fit.ci(2,:),'--r')
plot(limZ,fit.beta(2,:)-fit.ci(2,:),'--r')
ylim([0 70])
xlim([-700 -400])
set(gca,'XDir','reverse')
set(gca,'YScale','log')
xlabel('position [um], top \leftarrow \rightarrow bottom')
ylabel('{\itT}_2 [ms]')

