clear
clc
close all

%%

filename = 'Double_Gly_15mMGdH2O_T1IRBURP_29Sep2015';
filedir = 'C:\users\jnking01\desktop\';

fileloc = strcat(filedir,filename,'.tnt');

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tEcho = 700; %us
echoVector = (tEcho:tEcho:nEchoes*tEcho);

nEchoes = 64;
nPts = 76;
nPtsBlank = 2;
nT1Pts = 21;
T1min = 0.1; %ms
T1max = 40; %ms

linORlog = 1; % 0 for linearly space and 1 for log spaced

if linORlog == 0
    T1vector = linspace((T1min),(T1max),nT1Pts); % Linspace T1 points
else
    T1vector = logspace(log10(T1min),log10(T1max),nT1Pts); % Logspace T1sat
end

data = reshape(spec2',nPts,nEchoes,nT1Pts);
data = data(1:(nPts-nPtsBlank),:,:);
dataInt = sum(sum(data,1),2);
dataInt = reshape(dataInt,1,nT1Pts);
dataIntRe = real(dataInt);
dataIntIm = imag(dataInt);

%% Make 2D data set for T1IRT2 ILT

data2d = sum(real(data),1);
data2d = reshape(data2d,nEchoes, nT1Pts);
data2d = data2d';

surf(data2d); shading flat

save(strcat(filedir,filename,'.dat'), 'data2d', '-ascii')

%% 1D Fits

%T1
cftool(T1vector, dataIntRe./max(dataIntRe))

%T2
cftool(echoVector, data2d(:,end)./max(data2d))
