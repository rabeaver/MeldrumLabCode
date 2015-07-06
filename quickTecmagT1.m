clear
clc
close all

%%
filename = 'SBR_T1Sat_64ech_01Dec2014_15nD2.tnt';
filedir = 'C:\Users\tkmeldrum\Desktop\';

fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
% enter your values below
tEcho = 150; %us
nEchoes = 64;
nPts = 69;
nPtsBlank = 5;
nT1Pts = 15;
T1min = 0.390; %ms
T1max = 50; %ms
T1guess = 10 ; %ms 
T1lin = 0;
T1log = 1;
% stop entering your values

if T1lin == 1
    T1vector = linspace(T1min,T1max,nT1Pts);
elseif T1log == 1
    T1vector = linspace(log10(T1min),log10(T1max),nT1Pts);
    T1vector = 10.^(T1vector);
end
    
echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;    
    
data = reshape(spec2',nPts,nEchoes,nT1Pts);
data = data(1:(nPts-nPtsBlank),:,:);
dataInt = sum(sum(data,1),2);
dataInt = reshape(dataInt,1,nT1Pts);
dataIntRe = real(dataInt);
dataIntIm = imag(dataInt);
dataIntAb = abs(dataInt);

guesses = [0, max(dataIntRe), T1guess];
beta = nlinfit(T1vector,dataIntRe,@T1_recovery,guesses);
ypred = T1_recovery(beta,T1vector);

figure(1)
hold on
scatter(T1vector,dataIntRe);
plot(T1vector,ypred);


