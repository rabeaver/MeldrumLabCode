%% Incorporate all T2 values for solution, bead pack.
clear
clc
close all

tE = [50,60,70,80,90,100]*1e-6; % time (s) for different echo times
A1_bp = [0.3028	0.2981	0.2923	0.2813	0.3368	0.3153]; %A1, bead pack
A2_bp = [0.6655	0.6791	0.6874	0.6959	0.6455	0.6802]; %A2, bead pack
A_s = [0.9451	0.9498	0.9512	0.9478	0.9507	0.9552]; %A, solution
T2_1_bp = [2.8292	2.7266	2.5574	2.5117	2.8218	2.5212]; %T2, 1 bead pack (ms)
T2_2_bp = [10.3600	9.4823	8.6264	7.8357	7.3409	6.5706]; %T2, 2 bead pack (ms)
T2_s = [20.8898	17.4056	14.7492	12.3224	10.7296	9.0404]; %T2, solution (ms)
R2_1_bp = 1./T2_1_bp; %R2, 1 bead pack (ms-1)
R2_2_bp = 1./T2_2_bp; %R2, 2 bead pack (ms-1)
R2_s = 1./T2_s; %R2, solution (ms-1)

R2_eff_bp = A1_bp.*R2_1_bp + A2_bp.*R2_2_bp;
T2_eff_bp = 1./R2_eff_bp;

figure
subplot(2,2,1)
hold on
scatter(tE,A1_bp)
scatter(tE,A2_bp)
title('Amp, bead pack')
xlabel('echo time [s]')
ylim([0 1.1])
subplot(2,2,2)
hold on
scatter(tE,T2_1_bp)
scatter(tE,T2_2_bp)
title('T_2, bead pack')
xlabel('echo time [s]')
ylabel('T_2 [ms]')
ylim([0 30])
subplot(2,2,3)
scatter(tE,A_s)
title('Amp, solution')
xlabel('echo time [s]')
ylim([0 1.1])
subplot(2,2,4)
scatter(tE,T2_s)
title('T_2, solution')
xlabel('echo time [s]')
ylabel('T_2 [ms]')
ylim([0 30])

figure
scatter(tE,T2_eff_bp)

%% BIEXP FIT
clear 
close all
clc

tEref = 50; %us

datadir = '/Users/tyler/Desktop/Bead Pack/1500uMGd_BEADPACK_100u/1/';
datafile = 'data2.csv';
omitEchoes = 2;

parfilestem = strcat(datadir,'acqu');


% Load data and params
data = load(strcat(datadir,datafile));
params.acqTime = readpar_Kea(strcat(parfilestem,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(parfilestem,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(parfilestem,'.par'),'nrScans');
params.rxPhase = readpar_Kea(strcat(parfilestem,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(parfilestem,'.par'),'rxGain');
params.nrPts = readpar_Kea(strcat(parfilestem,'.par'),'nrPnts');
params.repTime = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.b1Freq = readpar_Kea(strcat(parfilestem,'.par'),'b1Freq');
params.nrEchoes = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoes');
params.echoTime = readpar_Kea(strcat(parfilestem,'.par'),'echoTime');

useEchoes = round(params.nrEchoes*(tEref/params.echoTime));

echoVector = params.echoTime:params.echoTime:params.nrEchoes*params.echoTime;
echoVector = echoVector'*1e-3; %switch to ms

% reshape data for formatting
data = reshape(data,(params.nrEchoes*2),params.nrPts);
dataRe = data(1:params.nrEchoes,:);
dataIm = data((params.nrEchoes+1):(params.nrEchoes*2),:);
dataCp = complex(dataRe,dataIm);
dataInt = sum(dataCp,2);
dataInt = dataInt./max(dataInt);

% plot(echoVector,real(dataInt))

% Fit to exp decay
A_guess = 0.9;
t2_guess = 3;
A_guess2 = 0.1;
t2_guess2 = 7; %ms

guesses = [A_guess;t2_guess;A_guess2;t2_guess2];

CI = 90; %desired confidence interval in percent

for i = 1:10
    [xfit,ypred,coeffs,coeffs_err,residuals] = bidecay_t2fit_simple(echoVector(omitEchoes+1:useEchoes),abs(dataInt(omitEchoes+1:useEchoes)),guesses,CI);
    guesses = coeffs;
end

figure
subplot(2,1,1)
hold on
plot(echoVector(omitEchoes+1:useEchoes),abs(dataInt(1+omitEchoes:useEchoes)))
plot(xfit,ypred)
plot(echoVector(omitEchoes+1:useEchoes), coeffs(1)*exp(-echoVector(omitEchoes+1:useEchoes)./coeffs(2)), '-r')
plot(echoVector(omitEchoes+1:useEchoes), coeffs(3)*exp(-echoVector(omitEchoes+1:useEchoes)./coeffs(4)), '-b')
subplot(2,1,2)
plot(echoVector(omitEchoes+1:useEchoes),residuals)

%% MONOEXP FIT
clear 
close all
clc

tEref = 50; %us

datadir = '/Users/tyler/Desktop/Bead Pack/1500uMGd_SOLUTION_100u/1/';
datafile = 'data2.csv';
omitEchoes = 2;
parfilestem = strcat(datadir,'acqu');

% Load data and params
data = load(strcat(datadir,datafile));
params.acqTime = readpar_Kea(strcat(parfilestem,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(parfilestem,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(parfilestem,'.par'),'nrScans');
params.rxPhase = readpar_Kea(strcat(parfilestem,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(parfilestem,'.par'),'rxGain');
params.nrPts = readpar_Kea(strcat(parfilestem,'.par'),'nrPnts');
params.repTime = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.b1Freq = readpar_Kea(strcat(parfilestem,'.par'),'b1Freq');
params.nrEchoes = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoes');
params.echoTime = readpar_Kea(strcat(parfilestem,'.par'),'echoTime');

useEchoes = round(params.nrEchoes*(tEref/params.echoTime));

echoVector = params.echoTime:params.echoTime:params.nrEchoes*params.echoTime;
echoVector = echoVector'*1e-3; %switch to ms

% reshape data for formatting
data = reshape(data,(params.nrEchoes*2),params.nrPts);
dataRe = data(1:params.nrEchoes,:);
dataIm = data((params.nrEchoes+1):(params.nrEchoes*2),:);
dataCp = complex(dataRe,dataIm);
dataInt = sum(dataCp,2);
dataInt = dataInt./max(dataInt);

% plot(echoVector,real(dataInt))

% Fit to exp decay
A_guess = 0.5;
t2_guess = 5; %ms

guesses = [A_guess;t2_guess];

CI = 90; %desired confidence interval in percent
for i=1:10
    [xfit,ypred,coeffs,coeffs_err,residuals] = monodecay_t2fit_simple(echoVector(omitEchoes+1:useEchoes),abs(dataInt(omitEchoes+1:useEchoes)),guesses,CI);
    guesses=coeffs;
end

figure
subplot(2,1,1)
hold on
plot(echoVector(omitEchoes+1:useEchoes),abs(dataInt(omitEchoes+1:useEchoes)))
plot(xfit,ypred)
subplot(2,1,2)
plot(echoVector(omitEchoes+1:useEchoes),residuals)
