close all
clc
clear
cd('Z:\Data\BMF\MortarDiffusion') % directory file
%load Tecmag data
[ap,spec1,spec2] = readTecmag('2015.tnt'); %reads tecmag data

%User-specified parameters
nrEchoes = 32;
nr2DPnts = ap.td(2); % reads number points from the aq params
blankPnts = 5; %number of "zero" points at ends of acq windows, normally 5
deltaMin = 15e-6;
deltaMax = 70e-6;
deltaStep = 9;
Td = 1000e-6; %constant Td time in s



%%% Other NMR parameters, time axes
G = 281.4345; %field gradient in MHz m-1
gammaHz = 42.576; %MHz T-1
gammaRad = 267.522e6; %T-1 s-1
delta = linspace(deltaMin,deltaMax,deltaStep); %different delta times in s
DELTA = Td - 2.*delta; %different DELTA times in s
G_T = G/gammaHz; %gradient in T m-1

%reshape data--not necessary, but more intuitive
spec3 = reshape(spec2',length(spec1)/nrEchoes,nrEchoes,nr2DPnts);
spec3 = spec3(1:size(spec3,1)-blankPnts,:,:);

%integrate each echo then each echo train
sumData = sum(sum(spec3,1),2);
sumData = real(reshape(sumData,nr2DPnts,1));
% 
% figure(1)
% plot(sumData,'k')
 
% EQUATION
%
% ln(I/I0) = -gammaRad^2 *   G_T^2 *     delta^2 * (DELTA + (2/3)*delta) * D
%             s-2 T-2        T2 m-2      s2              s                 m2 s-1
%
% END EQUATION

y = log(sumData./sumData(1))';
x = -gammaRad^2*G_T^2.*delta.^2.*(DELTA + (2/3)*delta);
% x1 = -gammaRad^2*G_T^2;
% x2 = delta.^2.*(DELTA + (2/3)*delta);
% x = x1*x2;

% figure(2)
% plot(x,y)

Diff = x'\y' %D in m2 s-1

dEach = y'./x';

gammaRadError = 0.0000063e6; % error in gyromagnetic ratio from wikipedia
G_Terror = G *0.002; % error in field gradient, estimated here to be .2%

sD = ((-2*y./(gammaRad^3*G_T^2.*delta.^2.*(DELTA + (2/3)*delta))).^2 * gammaRadError^2 + (-2*y./(gammaRad^2*G_T^3.*delta.^2.*(DELTA + (2/3)*delta))).^2 * G_Terror^2).^0.5;

% other way of doing the error
% [p,S] = polyfit(tau,y,1);
% fit= p(1)*tau+p(2);
% plot(tau,fit,'-r');
% [fitci(1).ypred,fitci(1).delta] = polyval(p,tau,S) %this needs to be ad1usted so that the % error in the slopes (as compared to real values can be determined)

% sx = sum(x);% this does a least squares error fit of the slope of the lines
% sy = sum(y);
% sxx = sum(x.^2);
% sxy = sum(y.*x);
% syy = sum(y.^2);
% N = length(x);
% D = N*sxx-sx^2;

% a0(1) = (sxx*sy-sx*sxy)/D; 
% a1(1) = (N*sxy-sx*sy)/D; % this will give you a value for diffusion differing from the one calculated above

% stdErra(1).a0 = sqrt((sxx*(syy-a1(1)*sxy-a0(1)*sy)/((N-2)*D)));
% stdErra(1).a1 = sqrt(N/sxx)*stdErra(1).a0;% gives you a standard deviation for the whole slope of the line
%% 


%%%%% IF WE ARE LOOKING FOR G FROM D
D = 2.2e-9; %m2 s-1
x = -gammaRad^2.*delta.^2.*(DELTA + (2/3)*delta)*D;

G2 = x'\y';
G = sqrt(G2)*gammaHz;

