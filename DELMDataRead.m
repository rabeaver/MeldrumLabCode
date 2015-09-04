% DELM data proc
% 2 Sept 2015, TKM

clear
clc
close all

%% load data, paramaters
% datadir = '/Users/tyler/Desktop/TOSEND_Genpurp_DELM_15Aug/1/';
% datadir = '/Users/tyler/Desktop/TOSEND_Laromer_slide17_DELM_28Aug/2/';
datadir = '/Users/tyler/Desktop/TOSEND_Laromer_slide17_DELM_28Aug/1/';
datafile = 'dataRe.dat';
paramsfile = 'acqu.par';

params.acqTime = readpar_Kea(strcat(datadir,paramsfile),'acqTime');
params.bandwidth = readpar_Kea(strcat(datadir,paramsfile),'bandwidth');
params.nScans = readpar_Kea(strcat(datadir,paramsfile),'nrScans');
params.rxPhase = readpar_Kea(strcat(datadir,paramsfile),'rxPhase');
params.rxGain = readpar_Kea(strcat(datadir,paramsfile),'rxGain');
params.nPts = readpar_Kea(strcat(datadir,paramsfile),'nrPnts');
params.repTime = readpar_Kea(strcat(datadir,paramsfile),'repTime');
params.b1Freq = readpar_Kea(strcat(datadir,paramsfile),'b1Freq');
params.nEchoes = readpar_Kea(strcat(datadir,paramsfile),'nrEchoes');
params.tE = readpar_Kea(strcat(datadir,paramsfile),'echoTime');

timefile = 'data.dat';

data = load(strcat(datadir,datafile));
tauTimes = load(strcat(datadir,timefile));
tauTimes = tauTimes(:,1);
tauTimes2 = tauTimes.^2;

echoVector = (params.tE:params.tE:params.nEchoes*params.tE)*1e-6;

data2 = reshape(data,length(tauTimes),params.nPts,params.nEchoes); %data2 sums each complex point in an echo to produce one value for each echo
data2 = sum(data2,2);
data2 = reshape(data2,length(tauTimes),params.nEchoes);

data3 = sum(data2,2); %data3 sums up all the echoes for one tau point
data3 = data3./max(data3);

%% plot DELM vs tau (left), vs. tau^2 (right)
figure(1)
subplot(1,2,1)
plot(tauTimes,data3,'-o');
xlabel('tau [ms]')
subplot(1,2,2)
plot(tauTimes2,data3,'-o');
xlabel('tau^2 [ms^2]')

%% fit to line
% choose indices for start and end range of the line in tau^2
startind = 15;
endind = 27;
% genpurp: 21-34; laromer 1: 15-27l Laromer 2: 14-27


%fit
[p,S] = polyfit(tauTimes2(startind:endind),data3(startind:endind),1);
[t2Line, delta] = polyval(p,tauTimes2,S);

%add to plot
figure(1)
subplot(1,2,2)
hold on
plot(tauTimes2,t2Line,'-r')
xlim([min(tauTimes2) max(tauTimes2)]);
ylim([0 1])

%% line error calc
Sx = sum(tauTimes2(startind:endind));
Sy = sum(data3(startind:endind));
Sxx = sum(tauTimes2(startind:endind).^2);
Syy = sum(data3(startind:endind).^2);
Sxy = sum(tauTimes2(startind:endind).*data3(startind:endind));
N = endind - startind + 1;
D = N*Sxx-Sx^2;

slope = p(1)
int = p(2);
sig_int = ((Sxx*(Syy-p(1)*Sxy-p(2)*Sy))/((N-2)*D))^0.5;
sig_slope = (N/Sxx)^0.5*sig_int