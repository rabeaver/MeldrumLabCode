clc
close all
clear 

%% ILT from data
[aq,sp1,sp2,sp3] = readTecmag('C:\Users\benjamin\Documents\Data\Diffusion\glycerol_liquid_n32_19June2014.tnt');
sp1 = real(sp1);
sp1 = reshape(sp1,32,256);
sp1 = sp1(1:size(sp1,1)-5,1:256);
specSum = sum(sp1);
specSum = specSum(5:end);

alpha = 5e8;% 5e8 for DI water, looking at the other ones 3e8 for Jamestown and Armory both
omitpoints = 1;
% sumData = sum(data1(:,:,i));

lowLim = 0.001e-3; %min(echoVector)/10000; %
hiLim = 1000e-3; %max(echoVector)/10;
% lowLim = 0.1e-3; % for quickset hydraulic
% hiLim = 100e-3; % for quickset hyraulic
 %divides by 2 to make the ilt move faster
nrEchoes = 256;
echoTime = 150e-6;
echoAxis = linspace(echoTime,echoTime*nrEchoes,nrEchoes);
echoAxis = echoAxis(1:length(specSum));
nrILTSteps = length(echoAxis);

[sample.spectrum,sample.tau,sample.chisq,sample.compte] = upnnlsmooth1D(specSum',echoAxis',  lowLim, hiLim, alpha ,  -1,  nrILTSteps);

figure(2)
semilogx(sample.tau,sample.spectrum)
title('ILT plots')
hold on


%% Diffusion data
[ap,spec1,spec2,spec3] = readTecmag('C:\Users\benjamin\Documents\Data\Diffusion\Glycerol_2D11_n128_s128_01July2014.tnt');


%User-specified parameters
nrEchoes = 128;
wantedNrEchoes = nrEchoes;% this is used for trimming the echoes if you want less to get better data and fits
nr2DPnts = ap.td(2);
nr3DPnts = ap.td(3);
blankPnts = 5; %number of "zero" points at ends of acq windows
Td = 3000e-6;
nrScans = 256;
echoTime = 150e-6; % in seconds
% DELTA = 500e-6; % for fixed Delta
deltaStep = 11; % this is the NUMBER of delta steps (delta points in excel spreadsheet)

deltaMin = 15e-6;
deltaMax = 500e-6;

% other params
gammaHz = 42.576; %MHz T-1
gammaRad = 267.522e6; %T-1 s-1
G = 281.4345; %field gradient in MHz m-1
G_T = G/gammaHz; %gradient in T m-1

spec2 = spec2(1:deltaStep,:);



%%% Other NMR parameters, time axes
delta = linspace(deltaMin,deltaMax,deltaStep); %different delta times in s
DELTA = Td - 2.*delta; %different DELTA times in s
% Td = DELTA + 2.*delta; % for different Diffusion times

%reshape data--not necessary, but more intuitive
% data = reshape(spec2',length(spec1)/nrEchoes,nrEchoes,nr2DPnts);% for 2d experiments

% for j = 1:nr3DPnts;

spec1 = reshape(spec2',ap.td(1),deltaStep); % reshapes into points, delta points, experiments
data1 = reshape(spec1(:,:),length(spec2)/nrEchoes,nrEchoes,deltaStep); % reshapes so that each echo can be summed
data1 = data1(1:size(data1,1)-blankPnts,1:wantedNrEchoes,:);

% sp2 = reshape(spec3(2,:),ap.td(2),ap.td(1));
% data2 = reshape(sp2',length(spec3)/nrEchoes,nrEchoes,nr2DPnts);
% data2 = data2(1:size(data2,1)-blankPnts,:,:);

%integrate each echo then each echo train
sumData = sum(sum(data1,1),2);
sumData = real(reshape(sumData,deltaStep,1));
sumDataPoint1 = real(sum(data1,1)/nrScans);% this gives integrated values for the points 

% sumData2 = sum(sum(data1,1),2);
% sumData2 = real(reshape(sumData2,nr2DPnts,1));
% sumDataPoint2 = real(sum(data1,1)/nrScans);% this gives integrated values for the points 

%%% exp fit
% guesses = [600,0.015,200,0.005]; % for quickset mortar%
guesses = [12000,.02]; % these work well for lime mortars 
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');
% sumData2 = real(sum(spec3,1)/64);
echoAxis = linspace(echoTime,echoTime * wantedNrEchoes,wantedNrEchoes);
% for j = 1:nr3DPnts;

for i = 1:deltaStep%size(sumData2,3)
    tau(i).beta = [0,0];
    while tau(i).beta~= guesses
    [tau(i).beta,tau(i).resid,tau(i).J] = nlinfit(echoAxis,sumDataPoint1(1,:,i),@t2monofit_simple,guesses,fitopts);
    tau(i).pred = t2monofit_simple(tau(i).beta,echoAxis);
    guesses=tau(i).beta; %reallocates guess to new result
    end
    ci = nlparci(tau(i).beta,tau(i).resid,'jacobian',tau(i).J);
    error(i).ci = ci; % gives confidence intervals
    
%     figure(i)
%     hold on
%     plot(echoAxis,tau(i).pred)
%     plot(echoAxis,sumDataPoint1(1,:,i),'-r')
end

% figure(1)
% plot(sumData,'k')
resid = tau(1).resid;
pred = tau(1).pred;
dataSum = sumDataPoint1(1,:,1);
amp(1) = tau(1).beta(1);
for i=2:deltaStep
    dataSum = [dataSum;sumDataPoint1(1,:,i)];
    resid = [resid,tau(i).resid];
    pred = [pred;tau(i).pred];
    amp = [amp,tau(i).beta(1)];
end
ampBoth = amp;
% samp1 = sample(1,:);
% samp2 = sample(2,:);
% samp3 = sample(3,:);
% samp4 = sample(4,:);
% samp5 = sample(5,:);
% samp6 = sample(6,:);
% samp7 = sample(7,:);
% samp8 = sample(8,:);
% samp9 = sample(9,:);
% samp10 = sample(10,:);
% samp11 = sample(11,:);
% samp12 = sample(12,:);
% samp13 = sample(13,:);
% samp14 = sample(14,:);
% samp15 = sample(15,:);
% samp16 = sample(16,:);
% samp17 = sample(17,:);
% samp18 = sample(18,:);
% samp19 = sample(19,:);
% samp20= sample(20,:);
% samp21 = sample(21,:);
% samp22 = sample(22,:);
% samp23 = sample(23,:);
% samp24 = sample(24,:);
% samp25 = sample(25,:);
% samp26 = sample(26,:);

% sample(k).dataSum = dataSum(:);
% sample(k).resid = resid(:);
% sample(k).pred = pred(:);
% sample(k).amp = amp(:);
% sample(k).amp2 = amp2(:); 
% EQUATION
%
% ln(I/I0) = -gammaRad^2 *   G_T^2 *     delta^2 * (DELTA + (2/3)*delta) * D
%             s-2 T-2        T2 m-2      s2              s                 m2 s-1
%
% END EQUATION

y = log(ampBoth./ampBoth(1)); % for amplitudes from the fits

y2 = log(sumData./sumData(1))'; % for summed Data
x = -gammaRad^2*G_T^2.*delta.^2.*(DELTA + (2/3)*delta);
% x1 = -gammaRad^2*G_T^2;
% x2 = delta.^2.*(DELTA + (2/3)*delta);
% x = x1*x2;

% figure(2)
% plot(x,y)

diff.summed = real(x'\y2'); %D in m2 s-1
diff.amp = real(x'\y');
% diff(k).each = y'/x';
% 
% gammaRadError = 0.0000063e6; % error in gyromagnetic ratio from wikipedia
% G_Terror = G *0.002; % error in field gradient, estimated here to be .2%
% 
% sD = ((-2*y./(gammaRad^3*G_T^2.*delta.^2.*(DELTA + (2/3)*delta))).^2 * gammaRadError^2 + (-2*y./(gammaRad^2*G_T^3.*delta.^2.*(DELTA + (2/3)*delta))).^2 * G_Terror^2).^0.5;

% other way of doing the error
% [p,S] = polyfit(tau,y,1);
% fit= p(1)*tau+p(2);
% plot(tau,fit,'-r');
% [fitci(1).ypred,fitci(1).delta] = polyval(p,tau,S) %this needs to be ad1usted so that the % error in the slopes (as compared to real values can be determined)

sx = sum(x);% this does a least squares error fit of the slope of the lines
sy = sum(y);
sxx = sum(x.^2);
sxy = sum(y.*x);
syy = sum(y.^2);
N = length(x);
D = N*sxx-sx^2;

a0 = (sxx*sy-sx*sxy)/D; 
a1 = (N*sxy-sx*sy)/D; % this will give you a value for diffusion differing from the one calculated above
% 
stdErra.a0 = sqrt((sxx*(syy-a1*sxy-a0*sy)/((N-2)*D)));
stdError = sqrt(N/sxx)*stdErra.a0;% gives you a standard deviation for the whole slope of the line
% end



% %%%%% IF WE ARE LOOKING FOR G FROM D
% D = 2.2e-9; %m2 s-1
% x = -gammaRad^2.*delta.^2.*(DELTA + (2/3)*delta)*D;
% 
% G2 = x'\y';
% G = sqrt(G2)*gammaHz
diff.amp
diff.summed

