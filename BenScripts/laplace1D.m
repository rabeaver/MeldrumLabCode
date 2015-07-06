close all
clear
clc
cd('C:\Users\bmfortman\Documents\Data\MortarDiffusion')
%% load Tecmag data
[ap,spec1,spec2] = readTecmag('Sample1_2D26_n64_23June2014.tnt');

%% User-specified parameters
nrEchoes = 256;
nr2DPnts = ap.td(2);
blankPnts = 5; %number of "zero" points at ends of acq windows

%% reshape data--not necessary, but more intuitive
spec3 = reshape(spec2',length(spec1)/nrEchoes,nrEchoes,nr2DPnts);
spec3 = spec3(1:size(spec3,1)-blankPnts,:,:);

%% integrate each echo
sumData = sum(spec3(:,:,1));

alpha = 5e7;
omitpoints = 1;
echoVector= (1:nrEchoes) * 150e-6; % nrEchoes * the echo time (in us)
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
