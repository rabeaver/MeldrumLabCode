close all
clc
clear
dir1 = ('C:\Users\benjamin\Documents\MortarDiffusion');
dir2 = ('C:\Users\bmfortman\Documents\Data\MortarDiffusion');

cd(dir2);% 
%load Tecmag data
[ap,spec1,spec2,spec3] = readTecmag('Jamestown_2D51_n256_s512_25June2014.tnt');

%User-specified parameters
nrEchoes = 256;
nr2DPnts = ap.td(2);
nr3DPnts = ap.td(3);
blankPnts = 5; %number of "zero" points at ends of acq windows
deltaMin = 15e-6;
deltaMax = 415e-6;
deltaStep = 51;
% Td = 1000e-6; %constant Td time in s
G = 281.4345; %field gradient in MHz m-1
nrScans = 512;
echoTime = 150e-6; % in seconds
DELTA = 500e-6; % for fixed Delta


%%% Other NMR parameters, time axes
gammaHz = 42.576; %MHz T-1
gammaRad = 267.522e6; %T-1 s-1
delta = linspace(deltaMin,deltaMax,deltaStep); %different delta times in s
% DELTA = Td - 2.*delta; %different DELTA times in s
Td = DELTA + 2.*delta; % for different Diffusion times
G_T = G/gammaHz; %gradient in T m-1

%reshape data--not necessary, but more intuitive
% data = reshape(spec2',length(spec1)/nrEchoes,nrEchoes,nr2DPnts);% for 2d experiments

for j = 1:nr3DPnts;

sp1 = reshape(spec3(j,:),ap.td(2),ap.td(1)); % select the 3d experiment you want right after spec3
data1 = reshape(sp1',length(spec3)/nrEchoes,nrEchoes,nr2DPnts);
data1 = data1(1:size(data1,1)-blankPnts,:,:);

% sp2 = reshape(spec3(2,:),ap.td(2),ap.td(1));
% data2 = reshape(sp2',length(spec3)/nrEchoes,nrEchoes,nr2DPnts);
% data2 = data2(1:size(data2,1)-blankPnts,:,:);

%integrate each echo then each echo train
sumData = sum(sum(data1,1),2);
sumData = real(reshape(sumData,nr2DPnts,1));
sumDataPoint1 = real(sum(data1,1)/nrScans);% this gives integrated values for the points 

% sumData2 = sum(sum(data1,1),2);
% sumData2 = real(reshape(sumData2,nr2DPnts,1));
% sumDataPoint2 = real(sum(data1,1)/nrScans);% this gives integrated values for the points 

%%% exp fit
guesses = [5000,.20,890,.4];% [1,10,1,10]; % ; %these work well for sample 7
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');
% sumData2 = real(sum(spec3,1)/64);
echoAxis = linspace(echoTime,echoTime * nrEchoes,nrEchoes);
% for j = 1:nr3DPnts;

for i = 1:deltaStep%size(sumData2,3)
    tau(i).beta = [0,0,0,0];
    while tau(i).beta~= guesses
    [tau(i).beta,tau(i).resid,tau(i).J] = nlinfit(echoAxis,sumDataPoint1(1,:,i),@t2bifit_simple,guesses,fitopts);
    tau(i).pred = t2bifit_simple(tau(i).beta,echoAxis);
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
amp2(1) = tau(1).beta(3);
for i=2:deltaStep
    dataSum = [dataSum;sumDataPoint1(1,:,i)];
    resid = [resid,tau(i).resid];
    pred = [pred;tau(i).pred];
    amp = [amp,tau(i).beta(1)];
    amp2 = [amp2,tau(i).beta(3)];
end
ampBoth = amp + amp2;
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

sample(j).dataSum = dataSum(:);
sample(j).resid = resid(:);
sample(j).pred = pred(:);
sample(j).amp = amp(:);
sample(j).amp2 = amp2(:); 
% EQUATION
%
% ln(I/I0) = -gammaRad^2 *   G_T^2 *     delta^2 * (DELTA + (2/3)*delta) * D
%             s-2 T-2        T2 m-2      s2              s                 m2 s-1
%
% END EQUATION

y = log(ampBoth./ampBoth(1)); % for amplitudes from the fits

% y = log(sumData./sumData(1))'; % for summed Data
x = -gammaRad^2*G_T^2.*delta.^2.*(DELTA + (2/3)*delta);
% x1 = -gammaRad^2*G_T^2;
% x2 = delta.^2.*(DELTA + (2/3)*delta);
% x = x1*x2;

% figure(2)
% plot(x,y)

diff(j).total = x'\y' %D in m2 s-1

diff(j).each = y'./x';
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
% end



% %%%%% IF WE ARE LOOKING FOR G FROM D
% D = 2.2e-9; %m2 s-1
% x = -gammaRad^2.*delta.^2.*(DELTA + (2/3)*delta)*D;
% 
% G2 = x'\y';
% G = sqrt(G2)*gammaHz
clear('x','y','ampBoth','amp','amp2','sp1','data1','sumData','sumDataPoint1','tau','sample')
end

%% Diffusion versus t^0.5
timeAxis = Td(1:end).^0.5;

dDiff = dEach1+dEach2./2;

figure(5)
hold on
scatter((Td(2:end).^0.5),diff(1).each(2:end))
% scatter(Td(2:end),dEach1(2:end))
% sV = -(dEach - 2.2e-9 /


%% T2 Data from CPMG
[ap,d1,d2,d3] = readTecmag('Jamestown_n32_nE512_s512_25June2014.tnt'); % reading the file, assuming in the same dir as above

nrEchoes = 512;
nr2DPnts = ap.td(2);
nr3DPnts = ap.td(3);
blankPnts = 5; %number of "zero" points at ends of acq windows
nrScans = 512;
echoTime = 150e-6;

d1 = real(d1)/nrEchoes; % takes only the real data, for each individual scan
sp1 = reshape(d1',length(d1)/nrEchoes,nrEchoes); % reshape to makes summing easier
data = sp1(1:size(sp1,1)-blankPnts,:); % cuts out the blank points from the tecmag
dataSum = sum(data);% sums each echo
echoAxis = linspace(echoTime,echoTime * nrEchoes, nrEchoes);

% exponential fitting

guesses = [5200,0.027,1500,0.005];
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');

[fit.beta,fit.resid,fit.J] = nlinfit(echoAxis(4:end),dataSum(4:end),@t2bifit_simple,guesses,fitopts);
[fit.pred] = t2bifit_simple(fit.beta,echoAxis(4:end));
% close all
figure(1)
hold on
plot(echoAxis(4:end),dataSum(4:end))
plot(echoAxis(4:end),fit.pred)

%
5181.7*0.0272 % this was some initial stuff with rho, don't use this script for laplace transform processing
ans+1783*0.0052
ans/(5181.7+1783.1)
1/ans/(4.987e-8)
rho = ans;

%% laplace transform
alpha = 5e7;
omitpoints = 1;
% sumData = sum(data1(:,:,i));

lowLim = 1e-3; %min(echoVector)/10000; %
hiLim = 100e-3; %max(echoVector)/10;
nrILTSteps = length(echoAxis);

[sample.spectrum,sample.tau,sample.chisq,sample.compte] = upnnlsmooth1D(real(dataSum)',echoAxis',  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
figure(2)
semilogx(sample.tau*3*rho,sample.spectrum)
