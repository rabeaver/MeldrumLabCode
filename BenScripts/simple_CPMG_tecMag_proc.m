clc
close all
clear all

%% T2 Data from CPMG
cd('C:\Users\benjamin\Documents\Data\MortarDiffusion');
[ap,d1,d2,d3] = readTecmag('Sample5_n69_nE32_s128_28July2014.tnt'); % reading the file, assuming in the same dir as above

nrEchoes = 512;
nr2DPnts = ap.td(2);
nr3DPnts = ap.td(3);
blankPnts = 5; %number of "zero" points at ends of acq windows
nrScans = 512;
echoTime = 150e-6;

guesses = [0.8,0.046,0.2,0.006];% for biexp fit
% guesses = [1,0.023];% for monoexp fit you don't need to change anything else

d1 = real(d1)/nrEchoes; % takes only the real data, for each individual scan
sp1 = reshape(d1',length(d1)/nrEchoes,nrEchoes); % reshape to makes summing easier
data = sp1(1:size(sp1,1)-blankPnts,:); % cuts out the blank points from the tecmag
dataSum = sum(data);% sums each echo
echoAxis = linspace(echoTime,echoTime * nrEchoes, nrEchoes);


dataNorm = dataSum./(max(dataSum));% normalizing data
if length(guesses) == 2
    [fit.beta,fit.resid,fit.J] = nlinfit(echoAxis(1:length(dataNorm)),dataNorm,@t2monofit_simple,guesses);
    [fit.pred] = t2monofit_simple(fit.beta,echoAxis(1:length(dataNorm)));
    ci = nlparci(fit.beta,fit.resid,'jacobian',fit.J);
    [fit.cil] = t2monofit_simple(ci(:,1),echoAxis(1:length(dataNorm)));
    [fit.ciu] = t2monofit_simple(ci(:,2),echoAxis(1:length(dataNorm)));
elseif length(guesses) == 4
    [fit.beta,fit.resid,fit.J] = nlinfit(echoAxis(1:length(dataNorm)),dataNorm,@t2bifit_simple,guesses);
    [fit.pred] = t2bifit_simple(fit.beta,echoAxis(1:length(dataNorm)));
    ci = nlparci(fit.beta,fit.resid,'jacobian',fit.J);
    [fit.cil] = t2bifit_simple(ci(:,1),echoAxis(1:length(dataNorm)));
    [fit.ciu] = t2bifit_simple(ci(:,2),echoAxis(1:length(dataNorm)));
else
end
%plotting figure

figure(1)
hold on
plot(echoAxis,dataNorm)
plot(echoAxis,fit.pred)
plot(echoAxis,fit.ciu)
plot(echoAxis,fit.cil)

legend('Normalized data','predicted fit','upper confidence interval','lower confidence interval')

fit.beta