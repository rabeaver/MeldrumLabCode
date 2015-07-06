close all
clc
clear
cd('C:\Users\bmfortman\Documents\Data\MortarDiffusion')%('C:\Users\benjamin\Documents\MortarDiffusion')
%load Tecmag data
[ap,spec1,spec2] = readTecmag('Sample1_2D21_n256_s512_24June2014.tnt');

%User-specified parameters
nrEchoes = 256;
nr2DPnts = ap.td(2);
blankPnts = 5; %number of "zero" points at ends of acq windows
deltaMin = 15e-6;
deltaMax = 415e-6;
deltaStep = 21;
Td = 1000e-6; %constant Td time in s
nrScans = 512;
echoTime = 150e-3;
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
sumDataPoint = real(sum(spec3,1)/nrScans);
%% exp fit
guesses = [239, -1/-3.6,5661,-1/-0.2529];% [1,10,1,10]; % ; %these work well for sample 7
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');

echoAxis = linspace(echoTime,echoTime*nrEchoes,nrEchoes);


for i = 1:deltaStep%size(sumData2,3)
    % tau(i).beta = [0,0,0,0];
    % while tau(i).beta~= guesses
    [tau(i).beta,tau(i).resid,tau(i).J] = nlinfit(echoAxis,sumDataPoint(1,:,i),@t2bifit_simple,guesses,fitopts);
    tau(i).pred = t2bifit_simple(tau(i).beta,echoAxis);
    % guesses=tau(i).beta; %reallocates guess to new result
    % end
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
sample = sumDataPoint(1,:,1);
amp(1) = tau(1).beta(1);
amp2(1) = tau(1).beta(3);
for i=2:deltaStep
    sample = [sample;sumDataPoint(1,:,i)];
    resid = [resid,tau(i).resid];
    pred = [pred;tau(i).pred];
    amp = [amp,tau(i).beta(1)];
    amp2 = [amp2,tau(i).beta(3)];
end
ampBoth = amp + amp2;

%% equation fitting 
% EQUATION
%
% ln(I/I0) = -gammaRad^2 *   G_T^2 *     delta^2 * (DELTA + (2/3)*delta) * D
%             s-2 T-2        T2 m-2      s2              s                 m2 s-1
%
% END EQUATION
y = log(ampBoth./ampBoth(1)); % for amplitudes from the fits

% y = log(sumData./sumData(1))'; % for summed data
x = -gammaRad^2*G_T^2.*delta.^2.*(DELTA + (2/3)*delta);
% x1 = -gammaRad^2*G_T^2;
% x2 = delta.^2.*(DELTA + (2/3)*delta);
% x = x1*x2;

% figure(2)
% plot(x,y)

Diff = x'\y' %D in m2 s-1

dEach = y'./x';
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
G = sqrt(G2)*gammaHz

%% Diffusion versus t^0.5
Tdiff = [530
570
610
650
690
730
770
810
850
890
930
970
1010
1050
1090
1130
1170
1210
1250
1290
1330];

timeAxis = (Tdiff * 10e-5).^0.5;

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
