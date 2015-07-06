clc
close all
clear all
%% Reading in data Jamestown
t2Data = ('Jamestown_n32_nE512_s128_1July2014.tnt');
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

JdataNorm = dataSum./max(dataSum);

diData = t2bifit_simple([0.8741,0.0565,0.0817,0.0057],echoAxis);
guesses = [0.25,0.4832,0.054,0.3,0.02]%,0.15,0.004];
lowerBounds = [0,0,0,0,0];
upperBounds = [1,inf,inf,inf,inf];
% fit.beta = [0,0,0,0,0,0,0];
% while fit.beta ~= guesses 
[lsqfit.beta,lsqfit.resid,lsqfit.J] =lsqcurvefit(@bifitH20scale,guesses,echoAxis,JdataNorm,lowerBounds,upperBounds);
guesses = lsqfit.beta(2:end);
waterPer = lsqfit.beta(1);
% end


Jdif= JdataNorm - waterPer*diData;

[fit.beta,fit.resid,fit.J] = nlinfit(echoAxis,Jdif,@t2bifit_simple,guesses);
fit.pred = t2bifit_simple(fit.beta,echoAxis);
ci = nlparci(fit.beta,fit.resid,'jacobian',fit.J);

figure(1)
plot(fit.pred)
hold on
plot(Jdif)
plot(diData)
legend('fit','sampledata','DiData')


%% reading in data Armory

t2Data = ('Armory_n32_nE512_s128_1July2014.tnt');
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

AdataNorm = dataSum./max(dataSum);

diData = t2bifit_simple([0.8741,0.0565,0.0817,0.0057],echoAxis);
guesses = [0.25,0.4832,0.054,0.3,0.02]%,0.15,0.004];
lowerBounds = [0,0,0,0,0];
upperBounds = [1,inf,inf,inf,inf];
% fit.beta = [0,0,0,0,0,0,0];
% while fit.beta ~= guesses 
[lsqfit.beta,lsqfit.resid,lsqfit.J] =lsqcurvefit(@bifitH20Scale,guesses,echoAxis,AdataNorm,lowerBounds,upperBounds);
guesses = lsqfit.beta(2:end);
waterPer = lsqfit.beta(1);
% end


Adif = AdataNorm - waterPer*diData;

[fit.beta,fit.resid,fit.J] = nlinfit(echoAxis,Adif,@t2bifit_simple,guesses);
fit.pred = t2bifit_simple(fit.beta,echoAxis);
ci = nlparci(fit.beta,fit.resid,'jacobian',fit.J);

figure(1)
plot(fit.pred)
hold on
plot(Adif)
plot(diData)
legend('fit','sampledata','DiData')

%%
alpha = 5e8;% 5e8 for DI water, looking at the other ones 3e8 for Jamestown and Armory both
omitpoints = 1;
% sumData = sum(data1(:,:,i));

lowLim = 1e-3; %min(echoVector)/10000; %
hiLim = 1000e-3; %max(echoVector)/10;
% lowLim = 0.1e-3; % for quickset hydraulic
% hiLim = 100e-3; % for quickset hyraulic
nrILTSteps = length(echoAxis); %divides by 2 to make the ilt move faster

[sample.spectrum,sample.tau,sample.chisq,sample.compte] = upnnlsmooth1D(JdataNorm',echoAxis',  lowLim, hiLim, alpha ,  -1,  nrILTSteps);

figure(2)
semilogx(sample.tau,sample.spectrum)
hold on

