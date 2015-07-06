close all
clc
clear
cd('/Users/jaredking/Documents/Research Files and Data/Paint/Diffusion')
%load Tecmag data
[ap,spec1,spec2] = readTecmag('Isopropanol_Diffusion.tnt');

%User-specified parameters
nrEchoes = 32;
nr2DPnts = ap.td(2);
blankPnts = 5; %number of "zero" points at ends of acq windows
deltaMin = 10e-6;
deltaMax = 70e-6;
deltaStep = 7;
Td = 150e-6; %constant Td time in s
G = 281.4345; %field gradient in MHz m-1


%%% Other NMR parameters, time axes
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
% 
% a0(1) = (sxx*sy-sx*sxy)/D; 
% a1(1) = (N*sxy-sx*sy)/D; % this will give you a value for diffusion differing from the one calculated above
% 
% stdErra(1).a0 = sqrt((sxx*(syy-a1(1)*sxy-a0(1)*sy)/((N-2)*D)));
% stdErra(1).a1 = sqrt(N/sxx)*stdErra(1).a0;% gives you a standard deviation for the whole slope of the line
%% 


%%%%% IF WE ARE LOOKING FOR G FROM D
D = 2.2e-9; %m2 s-1
x = -gammaRad^2.*delta.^2.*(DELTA + (2/3)*delta)*D;

G2 = x'\y';
G = sqrt(G2)*gammaHz;

%% Diffusion versus t^0.5
d1 = [9.000
12.000
15.000
18.000
21.000
24.000
27.000
30.000
33.000
36.000
39.000
42.000
45.000
48.000
51.000
54.000
57.000
60.000
63.000
66.000
69.000
72.000
75.000
78.000
81.000
84.000];

timeAxis = (d1 * 10e-5).^0.5;

plot(timeAxis(2:end),dEach(2:end))
% sV = -(dEach - 2.2e-9 /

%% laplace transform
clear
clc
close all
% try
%     matlabpool
% catch
% end
%% Load data
filename = 'T2data-1975.dat';
parname = 'acqu';
alpha = 5e7;
omitpoints = 1;
echoVector= (1:nrEchoes) * 150e-6;
for i = 1:deltaStep;
sumData = sum(spec3(:,:,i));
% sumData= sumData(:,1)

% data = load(filename);
% echoVector = data(omitpoints+1:end,1);
% real_data =  real(spec1);
% imag_data =  imag(spec1);

% real_data = load(strcat(filename,'-decaysRe.dat'));
% imag_data = load(strcat(filename,'-decaysIm.dat'));
% 
% contrast_data = load(strcat(filename,'.dat'));
% position = contrast_data(:,1);

% cplx_data = complex(real_data,imag_data);
% phased_data = autophase(cplx_data,0.1);
% 
% params.acqTime = readpar_Kea(strcat(parname,'.par'),'acqTime');
% params.bandwidth = readpar_Kea(strcat(parname,'.par'),'bandwidth');
% params.nrScans = readpar_Kea(strcat(parname,'.par'),'nrScansT2');
% params.rxPhase = readpar_Kea(strcat(parname,'.par'),'rxPhase');
% params.rxGain = readpar_Kea(strcat(parname,'.par'),'rxGain');
% params.nrPts = readpar_Kea(strcat(parname,'.par'),'nrPnts');
% params.repTime = readpar_Kea(strcat(parname,'.par'),'repTimeT2');
% params.b1Freq = readpar_Kea(strcat(parname,'.par'),'b1Freq');
% params.nrEchoes = readpar_Kea(strcat(parname,'.par'),'nrEchoesT2');
% params.echoTime = readpar_Kea(strcat(parname,'.par'),'echoTime');

% echoVector = params.echoTime:params.echoTime:params.echoTime*params.nrEchoes;
lowLim = 0.1e-3; %min(echoVector)/10000; %
hiLim = 10000e-3; %max(echoVector)/10;
nrILTSteps = length(echoVector);

[sample(i).spectrum,sample(i).tau,sample(i).chisq,sample(i).compte] = upnnlsmooth1D(real(sumData)',echoVector',  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
figure(i)
semilogx(sample(i).tau,sample(i).spectrum)
end
