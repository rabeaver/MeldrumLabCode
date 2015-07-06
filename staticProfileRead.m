clear
close all
clc
%% User parameters

% Input data location and filestem, number of profiles, depths per profile,
% and echoes per depth

cd('C:\Users\TMeldrum\Dropbox\Data from MobileMOUSE\CdYellow\Profiles_Day2ofDrying\threedayexpt\')
filestem = 'staticprof_CdS';

params.numberProfiles = readpar_Kea(strcat(filestem,'.par'),'nrExp');
params.finalDepth = readpar_Kea(strcat(filestem,'.par'),'finalDepth');
params.initDepth = readpar_Kea(strcat(filestem,'.par'),'initDepth');
params.stepSize = readpar_Kea(strcat(filestem,'.par'),'stepSize');
params.numberDepths = 1 + (params.finalDepth - params.initDepth)/params.stepSize;
params.numberEchoes = readpar_Kea(strcat(filestem,'.par'),'nrEchoes');
params.rxPhase = readpar_Kea(strcat(filestem,'.par'),'rxPhase');
params.bandwidth = readpar_Kea(strcat(filestem,'.par'),'bandwidth');
params.pulseLength = readpar_Kea(strcat(filestem,'.par'),'pulseLength');
params.b1Freq = readpar_Kea(strcat(filestem,'.par'),'b1Freq');
params.numberScans = readpar_Kea(strcat(filestem,'.par'),'nrScans');
params.repTime = readpar_Kea(strcat(filestem,'.par'),'repTime');
params.echoTime = readpar_Kea(strcat(filestem,'.par'),'echoTime0');
% 
% params.numberProfiles = 353;
% params.finalDepth = 4000;
% params.initDepth = 4000;
% params.stepSize = -1;
params.numberDepths = 1; % + (params.finalDepth - params.initDepth)/params.stepSize;
% params.numberEchoes = 6;
% params.rxPhase = 81;
% params.bandwidth = 1000;
% params.pulseLength = 25;
% params.b1Freq = 13.76;
% params.numberScans = 150;
% params.repTime = 200;
% params.echoTime = 115;
%% Load data

% The two lines with echoTime grab the time points for the echoes in ms
% from file 1.
% The realData takes the real data from the n profiles and stores them in
% the large array realData. Same for imaginary data. Check that the files
% end with "xx-decaysRe.dat" and are tab-delimited.

exptime = dlmread(strcat(filestem,'exptime.dat'));

if size(params.numberDepths) < 1;
    params.numberDepths = 1;
end

if params.numberProfiles > 1.56;
    depth = dlmread(strcat(filestem,'01.dat'));
    depth = depth(:,1);

    echoTime = dlmread(strcat(filestem,'01-decaysRe.dat'));
    echoTime = echoTime(:,1);

    realData = zeros(params.numberProfiles,params.numberDepths,params.numberEchoes);
%     normVal = zeros(params.numberProfiles,1);
    for i=1:1:params.numberProfiles
        temp = dlmread(strcat(filestem,num2str(i,'%02g'),'-decaysRe.dat'));
%         normVal(i) = max(temp(:,2:end));
        realData(i,:,:) = temp(:,2:end)'; %./normVal(i);
    end
    
    imagData = zeros(params.numberProfiles,params.numberDepths,params.numberEchoes);
    for i=1:1:params.numberProfiles
        temp = dlmread(strcat(filestem,num2str(i,'%02g'),'-decaysIm.dat'));
        imagData(i,:,:) = temp(:,2:end)'; %./normVal(i);
    end
else
    depth = dlmread(strcat(filestem,'.dat'));
    depth = depth(:,1);
    
    echoTime = dlmread(strcat(filestem,'-decaysRe.dat'));
    echoTime = echoTime(:,1);
    
    realData = dlmread(strcat(filestem,'-decaysRe.dat'));
    realData = realData(:,2:end);
    
    imagData = dlmread(strcat(filestem,'-decaysIm.dat'));
    imagData = imagData(:,2:end);
end
%% Autophase and save (speeds up for next time)
cplxData = complex(reshape(realData,params.numberProfiles,params.numberEchoes),reshape(imagData,params.numberProfiles,params.numberEchoes));
phasedData = autophase(cplxData,1)';
save('allPhasedData.mat','phasedData','cplxData','echoTime','exptime','params')
%% Load saved data (from previous cell)
load('C:\Users\TMeldrum\Dropbox\Data from MobileMOUSE\CdYellow\Profiles_Day2ofDrying\threedayexpt\allPhasedData.mat')
rawArea = sum(real(phasedData),1);
% for i = 1:length(rawArea)
%     phasedData(:,i) = phasedData(:,i)./rawArea(i);
% end
%%
blockSize = 1;
blockData = zeros(length(echoTime),length(exptime)/blockSize);
realData = real(phasedData);
for i = blockSize:blockSize:length(exptime)
    blockData(:,i/blockSize) = sum(realData(:,i-blockSize+1:i),2);
end
newExpTime = exptime(blockSize:blockSize:end);
%% Do ILT

alpha = 4e9;
lowLim = 0.05;
hiLim = 5000;
nrILTSteps = 200;
clear('spectrum','tau','chisq')
spectrum = zeros(nrILTSteps,length(exptime)/blockSize);
tau = spectrum;
chisq = zeros(length(exptime)/blockSize,1);

matlabpool(4)
tic
parfor i = 1:size(blockData,2)
%     tic
    [spectrum(:,i),tau(:,i),chisq(i),~     ] = upnnlsmooth1D(blockData(:,i),echoTime,  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
    disp(strcat('No. ',num2str(i)));
%     toc
end
totalTime = toc;
disp(strcat('Time per loop: ',num2str(totalTime/(params.numberProfiles/blockSize)),' sec.'));
matlabpool close

ILTArea = sum(spectrum,1);
% rawArea = sum(real(phasedData),1);
% save('group7.mat','spectrum','tau')

% figure
% semilogx(tau,spectrum)
%% Plot
figure(1)
subplot(2,2,[1 2])
hold on
surf(newExpTime'/3600,log10(tau),spectrum);
colormap(hot)
 shading interp
 set(gca,'YDir','Reverse')
 xlabel('time (hours)')
 ylabel('log(T_2/ms)')
 zlabel('intensity')
 ylim([log10(lowLim) log10(hiLim)])
 xlim([0 max(newExpTime)/3600])
 line([0;max(newExpTime)],[1;1],'Color','w')
 line([0;max(newExpTime)],[0;0],'Color','w')
 line([0;max(newExpTime)],[2;2],'Color','w')
 line([0;max(newExpTime)],[-1;-1],'Color','w')
 line([0;max(newExpTime)],[3;3],'Color','w')
 subplot(2,2,3)
 plot(newExpTime/3600,ILTArea,'-b')
 ylabel('Total area of ILT curve')
 xlim([0 max(newExpTime)/3600])
 xlabel('time (hours)')
 subplot(2,2,4)
 plot(exptime/3600,rawArea,'-k')
 xlim([0 max(exptime)/3600])
 xlabel('time (hours)')
 ylabel('Raw area of decay curve')