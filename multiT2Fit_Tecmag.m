% multiT2Fit_Tecmag.m
% 15 July 2015, TKM
% Fits multiple CPMG measurements (from a single 2D Tecmag file) to
% biexponential fits, outputs the corresponding T2 and amplitude values

clear
clc
close all

%% Load data
if ispc == 1;
    filedir = 'C:\Users\tkmeldrum\Desktop\LaromerUA9033\'; 
    filename = 'MilkT2D_2.tnt';
elseif ismac == 1;
    filedir = '/Users/tyler/Desktop/LaromerUA9033/';
    filename = 'LaromerUA9033_6%wtPI_slide12_CPMG_01_airdry_09Jul.tnt';
end

[aq,spec,spec2] = readTecmag4d(strcat(filedir,filename));

n2D = aq.td(2);
nPts = aq.ctd;
nEchoes = aq.td(1)/nPts;
tE = 150e-6; %s echo time
echoVector = linspace(tE, tE*nEchoes, nEchoes);
omitPts = 5; % zero points to clear at the end of each acq period

spec3 = reshape(spec2,n2D,nPts,nEchoes);
spec3 = spec3(:,1:end-omitPts,:); %removes zero points

data = sum(spec3,2); %adds (integrates) each echo...
data = reshape(data,n2D,nEchoes); %... and reshapes to normal
data = data./(max(max(data))); %normalizes to 1

%% Fit to monoexp

monoguesses = [0.5;
               50e-6]; % guesses (A; T2)
CI = 90;         

for i = 1:n2D
    [xfit,ypred(i,:),beta(:,i),beta_err(:,i),resid(:,i)] = monodecay_t2fit_simple(echoVector,abs(data(i,:)),monoguesses,CI);
    monoguesses = beta(:,i); %set guesses for next iteration to previous output
end