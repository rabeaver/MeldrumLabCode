clear
clc
close all

%%
% This is the file from the Kea spectrometer that you need to load. It
% reads in two columns: the first has the echo times in microseconds (us),
% and the second is the signal intensity at that point. You should select
% alpha manually, uaually between 1e6 and 1e10.

filedir = '/Users/tyler/Dropbox/Data/NGA/NGA_5June2015/Atlas9pARLO_1024/1/';

cd(filedir)
filename = 'data.csv';
parname = 'acqu';
alpha = 1e9;
omitpoints = 0;

data = load(filename);
echoVector = data(omitpoints+1:end,1)/1e6;
realData = data(omitpoints+1:end,2)./max(data(omitpoints+1:end,2));

lowLim = 10^-4; %s
hiLim = 10^-2;  %s
nrILTSteps = min(128,length(echoVector));

scatter(echoVector,realData)

kernel1 = 'exp(-h/T)';
[spectrum,tau] = upnnlsmooth1D(realData,echoVector,  lowLim, hiLim, alpha ,  -1,  nrILTSteps,kernel1);

spectrum = spectrum./sum(spectrum);

figure
hold on
yyaxis left
semilogx(tau,spectrum)
xlabel('T_2 time [s]')
ylabel('intensity [arb]')

yyaxis right
plot(tau,cumsum(spectrum))
ylim([0 1.05])
ylabel('integrated intensity')
set(gca,'XScale','log')


ILTout = [tau',spectrum'];

% save('ILTout.dat','ILTout','-ascii');