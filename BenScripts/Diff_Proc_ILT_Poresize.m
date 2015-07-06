close all
clc
clear
%% 
% compDir = ('C:\Users\bmfortman\Documents\Data');% lab Compy
compDir = ('C:\Users\benjamin\Documents\Data');% personal laptop
dir1 = strcat(compDir,'\MortarDiffusion');

diffData = ('Sample1_2D81_n128_s512_24July2014.tnt');
t2Data = ('Armory_n69_nE512_s128_21July2014.tnt');
saveData = 0; % 1if you want it to save ILT data, 0 if not
saveName = ('Sample7'); %the first part of the name it will be saving it as

cd(dir1);% 
%load Tecmag data
[ap,spec1,spec2,spec3] = readTecmag(diffData);

%User-specified parameters
nrExp = 9;
nrEchoes = 128;
wantedNrEchoes = nrEchoes;% this is used for trimming the echoes if you want less to get better data and fits
nr2DPnts = ap.td(2);
nr3DPnts = ap.td(3);
blankPnts = 5; %number of "zero" points at ends of acq windows
Td = linspace(1000e-6,1800e-6,nrExp); %Td time in s for each experiment
nrScans = 256;
echoTime = 150e-6; % in seconds
% DELTA = 500e-6; % for fixed Delta
deltaStep = 9; % this is the NUMBER of delta steps (delta points in excel spreadsheet)


% other params
gammaHz = 42.576; %MHz T-1
gammaRad = 267.522e6; %T-1 s-1
G = 281.4345; %field gradient in MHz m-1
G_T = G/gammaHz; %gradient in T m-1

spec2 = spec2(1:nrExp*deltaStep,:);

for k = 1:nrExp;

deltaMin = Td(k) * 0.1;
deltaMax = Td(k) * 0.4;

%%% Other NMR parameters, time axes

delta = linspace(deltaMin,deltaMax,deltaStep); %different delta times in s
DELTA = Td(k) - 2.*delta; %different DELTA times in s
% Td = DELTA + 2.*delta; % for different Diffusion times

%reshape data--not necessary, but more intuitive
% data = reshape(spec2',length(spec1)/nrEchoes,nrEchoes,nr2DPnts);% for 2d experiments

% for j = 1:nr3DPnts;

sp1 = reshape(spec2',ap.td(1),deltaStep,nrExp); % reshapes into points, delta points, experiments
data1 = reshape(sp1(:,:,k),length(spec2)/nrEchoes,nrEchoes,deltaStep); % reshapes so that each echo can be summed
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
guesses = [5000,.20,890,.4]; % these work well for lime mortars 
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');
% sumData2 = real(sum(spec3,1)/64);
echoAxis = linspace(echoTime,echoTime * wantedNrEchoes,wantedNrEchoes);
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

sample(k).dataSum = dataSum(:);
sample(k).resid = resid(:);
sample(k).pred = pred(:);
sample(k).amp = amp(:);
sample(k).amp2 = amp2(:); 
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

diff(k).summed = x'\y2'; %D in m2 s-1
diff(k).amp = x'\y';
%diff(k).each = y'./x';
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

dEach = diff(1).summed;
for l= 2:nrExp;
    
dEach = [dEach;diff(l).summed]; %converts from structure to a matrix
end
pointTrunc = 2;
[a,b] = polyfit(timeAxis(pointTrunc:end)',dEach(pointTrunc:end),1);
sv = a(1)/(-a(2)^1.5*(4/9/pi^0.5));
fitLine = a(1)*(Td(pointTrunc:end).^0.5) + a(2);

fitSlope= num2str(real(a(1)));% I got clever with my figures ;)
fitInt = num2str(real(a(2)));% oh yes, CLEVER.
fitLegend = strcat(saveName,'Y = ',fitSlope,'x + ',fitInt);% where the magic happens
figure(5)
hold on
scatter((Td(pointTrunc:end).^0.5),dEach(pointTrunc:end))
plot((Td(pointTrunc:end).^0.5),fitLine);
legend(fitLegend) %BAM it displays
xlabel ('square root time')
ylabel('diff value m^2/2')
% scatter(Td(2:end),dEach1(2:end))
% sV = -(dEach - 2.2e-9 /


%% T2 Data from CPMG
[ap,d1,d2,d3] = readTecmag(t2Data); % reading the file, assuming in the same dir as above

nrEchoes = 512;
nr2DPnts = ap.td(2);
nr3DPnts = ap.td(3);
blankPnts = 5; %number of "zero" points at ends of acq windows
nrScans = 128;
echoTime = 150e-6;

d1 = real(d1)/nrEchoes; % takes only the real data, for each individual scan
sp1 = reshape(d1',length(d1)/nrEchoes,nrEchoes); % reshape to makes summing easier
data = sp1(1:size(sp1,1)-blankPnts,:); % cuts out the blank points from the tecmag
dataSum = sum(data);% sums each echo
echoAxis = linspace(echoTime,echoTime * nrEchoes, nrEchoes);

echoAxisCon = echoAxis(1);
dataSumCon = dataSum(1);
for i = linspace(3,nrEchoes-1,length(echoAxis)/2-1)%cuts down nrPoints in the spectrum to allow them to work with the twoDlaplace
    echoAxisCon = [echoAxisCon; echoAxis(i)];
    dataSumCon = [dataSumCon; dataSum(i)];

end
dataNorm = dataSum./(max(dataSum));

% save('echoAxis.txt','echoAxisCon','-ascii'); % to save the data to be analyzed with the 2dLaplace program
% save('summedData.txt','dataSumCon','-ascii');
% TwoDLaplaceInverse
%% normalizing Data

% exponential fitting
guesses = [0.8,0.027,0.2,0.005];% for lime mortar
% guesses = [0.9,0.003,0.14,0.00093];% for hydraulic mortar

lowerBounds = [0,0,0,0];
upperBounds = [1,1,1,1];
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');

[lsqfit.beta,lsqfit.resid,lsqfit.J] =lsqcurvefit(@t2bifit_simple,guesses,echoAxis,dataNorm,lowerBounds,upperBounds);
guesses = lsqfit.beta(1:end); %takes results of lsq curve fit to use as guesses for nlinfit
% waterPer = lsqfit.beta(1); % water Percentage from mortar

% diData = t2bifit_simple([0.8741,0.0565,0.0817,0.0057],echoAxis);
% dataDif = dataNorm - waterPer*diData; % removes water weighted from data
% 
% dataDif = dataDif./max(dataDif); % renormalizes data after the subtraction occurs

dataDif = dataNorm; % when not removing the %water, just makes it the same as dataNorm for ease of use

[fit.beta,fit.resid,fit.J] = nlinfit(echoAxis(1:length(dataNorm)),dataDif,@t2bifit_simple,guesses,fitopts);
[fit.pred] = t2bifit_simple(fit.beta,echoAxis(1:length(dataDif)));
ci = nlparci(fit.beta,fit.resid,'jacobian',fit.J);
[fit.cil] = t2bifit_simple(ci(:,1),echoAxis(1:length(dataDif)));
[fit.ciu] = t2bifit_simple(ci(:,2),echoAxis(1:length(dataDif)));


% close all
figure(1)
hold on
plot(echoAxis(1:length(dataDif)),dataDif)
plot(echoAxis(1:length(dataDif)),fit.pred)
plot(echoAxis(1:length(dataDif)),fit.ciu,'-r')
plot(echoAxis(1:length(dataDif)),fit.cil,'-k')
legend('data','pred','upperci','lowerci')

%creating box of t2 dist as errorbars
t21A = fit.beta(1)/(fit.beta(1)+fit.beta(3));%relative abundance of the first t2 based on amp
t21 = linspace(ci(2,1),ci(2,2),100);% creates a scale of t2 values, based on ci
t21AV = ones(length(t21),1).*t21A; % creates the vector that's the same size as t21
t22A = fit.beta(3)/(fit.beta(1)+fit.beta(3));% relative abun of the 2nd t2 based on amp
t22 = linspace(ci(4,1),ci(4,2),100);% creates a scale of t2 values based on ci
t22AV = ones(length(t21),1).*t22A;% creates the vector that's the same size as t21


 % uses harmonic mean and s-v ratio to find value for surface relaxivity
%
% 5181.7*0.0272
% ans+1783*0.0052
% ans/(5181.7+1783.1)
% 1/ans/(4.987e-8)
% rho = ans;

%% laplace transform
alpha = 5e8;% 5e8 for DI water, looking at the other ones 3e8 for Jamestown and Armory both
omitpoints = 1;
% sumData = sum(data1(:,:,i));

lowLim = 0.001e-3; %min(echoVector)/10000; %
hiLim = 10000e-3; %max(echoVector)/10;
% lowLim = 0.1e-3; % for quickset hydraulic
% hiLim = 100e-3; % for quickset hyraulic
nrILTSteps = length(echoAxis); %divides by 2 to make the ilt move faster

[sample.spectrum,sample.tau,sample.chisq,sample.compte] = upnnlsmooth1D(real(dataSum)',echoAxis',  lowLim, hiLim, alpha ,  -1,  nrILTSteps);

figure(2)
semilogx(sample.tau,sample.spectrum)
hold on

sumN = sum(sample.spectrum);% sums the y's
sumD = sum(sample.spectrum./sample.tau);% sums the weights over the x's

meanT2h = sumN/sumD; % gives the harmonic mean of the distribution

rho = 1/meanT2h/sv; % gives surface relaxivity
% semilogx(sample.tau,sample.spectrum)
%% plotting ILT
intILT = cumsum(sample.spectrum);
intILT = intILT./max(intILT);
normILT = sample.spectrum./max(sample.spectrum);
if  t21A > t22A % multiplies the ILT by the higher of the two t2 abundances to give scaled values
    normILT = normILT * t21A; % multiplied by the relative t2 abundance to give values of the curves
else
    normILT = normILT * t22A; % multiplied by the relative t2 abundance to give values of the curves
end

figure(4)% plots ILT dist
semilogx(sample.tau,intILT,'-k')
hold on
semilogx(sample.tau,normILT,'-r')
plot(t21,t21AV,'-m','LineWidth',2)
plot(t22,t22AV,'-c','LineWidth',2)
stem(fit.beta(2),t21A,'-m')
stem(fit.beta(4),t22A,'-c')
xlabel('t2 times [s]')
ylabel('relative percent abundance')
legend('Integrated','Distribution','first T2 value from fit','2nd T2 Value from fit')


figure(3)% plots pore dist
xScalar = 1e6*3*rho*2;% the scaling factor that converts t2 times into pore diameter
semilogx(xScalar*sample.tau,intILT,'-k')
hold on
semilogx(xScalar*sample.tau,normILT,'-r')
plot(xScalar*t21,t21AV,'-m','LineWidth',2)
plot(xScalar*t22,t22AV,'-c','LineWidth',2)
stem(xScalar*fit.beta(2),t21A,'-m')
stem(xScalar*fit.beta(4),t22A,'-c')
xlabel('Pore Diameter [um]')
legend('Integrated','Distribution','first T2 value from fit','2nd T2 Value from fit')
set(gca,'XScale','log')

t21x = xScalar*t21;
t22x = xScalar*t22;
t21Ax = xScalar*fit.beta(2);
t22Ax = xScalar*fit.beta(4);

Xaxis = sample.tau*xScalar;
%% saving data
if saveData == 1;
save(strcat(saveName,'X.txt'),'Xaxis','-ascii');
save(strcat(saveName,'int.txt'),'intILT','-ascii');
save(strcat(saveName,'norm.txt'),'normILT','-ascii');
save(strcat(saveName,'rho.txt'),'rho','-ascii');
save(strcat(saveName,'t21AV.txt'),'t21AV','-ascii');
save(strcat(saveName,'t22AV.txt'),'t22AV','-ascii');
save(strcat(saveName,'t21A.txt'),'t21A','-ascii');
save(strcat(saveName,'t22A.txt'),'t22A','-ascii');
save(strcat(saveName,'t21x.txt'),'t21x','-ascii');
save(strcat(saveName,'t22x.txt'),'t22x','-ascii');
save(strcat(saveName,'t21Ax.txt'),'t21Ax','-ascii');
save(strcat(saveName,'t22Ax.txt'),'t22Ax','-ascii');

elseif saveData == 0;
    display('notsaved data = lost data')
else
end
