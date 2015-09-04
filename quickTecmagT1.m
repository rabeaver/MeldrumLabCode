clear
clc
close all

%%

filename = 'DoubleSample_15mMGdH2O_Glycerol_2DT1IR_BURP_LONGLONGLONG_30Aug2015.tnt';
filedir = 'C:\users\jnking01\desktop\messyprocfolder\';

fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tEcho = 700; %us

nEchoes = 64;
nPts = 76;
nPtsBlank = 2;
nT1Pts = 153;
T1min = 1.331; %ms
T1max = 31.681; %ms
T1guess = 0; %ms 

% T1vector = linspace((T1min),(T1max),nT1Pts); % Linspace T1sat
echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

T1vector = logspace(log10(T1min),log10(T1max),nT1Pts); % Logspace T1sat

data = reshape(spec2',nPts,nEchoes,nT1Pts);
data = data(1:(nPts-nPtsBlank),:,:);
dataInt = sum(sum(data,1),2);
dataInt = reshape(dataInt,1,nT1Pts);
dataIntRe = real(dataInt);
dataIntIm = imag(dataInt);


%% cftool
cftool(T1vector, -dataIntRe./min(dataIntRe))

%%
guesses = [1, max(dataIntRe), T1guess];
[beta, Resids, J] = nlinfit(T1vector,dataIntRe,@T1_recovery,guesses);
ypred = T1_recovery(beta,T1vector);
CI = nlparci(beta, Resids, 'jacobian', J);


figure(1)
hold on
scatter(T1vector,dataIntRe);
plot(T1vector,ypred);

%% Make 2D data set for T1SRT2

data2d = sum(real(data),1);
data2d = reshape(data2d,nEchoes, nT1Pts);
data2d = data2d';

surf(data2d); shading flat

save('DoubleSample_15mMGdH2O_Glycerol_2DT1IR_BURP_LONGLONGLONG_30Aug2015.dat', 'data2d', '-ascii')
