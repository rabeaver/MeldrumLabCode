clear
clc
close all
try
    matlabpool(4)
catch
end


%%
dataDir = '/Users/tyler/Dropbox/Data/Paint/CHARISMA Painting Study/Round 3/';
parfilestem = 'Spot4_profile';
omitpoints = 1;
interpFactor = 1;

cd(dataDir);
%% Load parameters
params.acqTime = readpar_Kea(strcat(parfilestem,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(parfilestem,'.par'),'bandwidth');
params.nrScansT1 = readpar_Kea(strcat(parfilestem,'.par'),'nrScans');
params.nrScansT2 = readpar_Kea(strcat(parfilestem,'.par'),'nrScansT2');
params.rxPhase = readpar_Kea(strcat(parfilestem,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(parfilestem,'.par'),'rxGain');
params.nrPtsT1 = readpar_Kea(strcat(parfilestem,'.par'),'nrPntsT1');
params.t1Max = readpar_Kea(strcat(parfilestem,'.par'),'tMax');
params.t1Est = readpar_Kea(strcat(parfilestem,'.par'),'t1Est');
params.repTimeT1 = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.repTimeT2 = readpar_Kea(strcat(parfilestem,'.par'),'repTimeT2');
params.b1Freq = readpar_Kea(strcat(parfilestem,'.par'),'b1Freq');
params.stepSize = readpar_Kea(strcat(parfilestem,'.par'),'stepSize');
params.finalDepth = readpar_Kea(strcat(parfilestem,'.par'),'finalDepth');
params.initDepth = readpar_Kea(strcat(parfilestem,'.par'),'initDepth');
params.nrEchoesT2 = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoesT2');
params.nrEchoes = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoes');
params.smartScan = readpar_Kea(strcat(parfilestem,'.par'),'smartScan');
%% Load depths

filename = 'Spot4_profile02';
% depthList = load('depths.dat');
depthList = load(strcat(filename,'.dat'));
depthList = depthList(:,1);
% depthList = 2500:-15:745;

% cd rough
% for n = 1:length(depthList)
%     temp = load(strcat('T2data_rough-',num2str(depthList(n)),'.dat'));
%     echoTime = temp(:,1);
%     realData(:,n) = temp(:,2);
%     imagData(:,n) = temp(:,3);
% end
% cd(dataDir)

realData = load(strcat(filename,'-decaysRe.dat'));
imagData = load(strcat(filename,'-decaysIm.dat'));
echoTime = realData(:,1);
realData = realData(:,2:end);
imagData = imagData(:,2:end);

cplxData = autophase(complex(realData,imagData),1);
echoTime = echoTime(omitpoints+1:end);
cplxData = cplxData(omitpoints+1:end,:);

realEndData = real(cplxData);

if interpFactor ~= 0
    echoTimeI = interp1(linspace(1,max(echoTime),length(echoTime)),echoTime,linspace(1,max(echoTime),length(echoTime)*interpFactor))';
    parfor j = 1:size(realEndData,2)
        realEndDataI(:,j) = interp1(linspace(1,max(echoTime),length(echoTime)),realEndData(:,j),linspace(1,max(echoTime),length(echoTime)*interpFactor));
    end
else
    echoTimeI = echoTime;
    realEndDataI = realEndData;
end
    

figure(1)
surf(depthList,echoTimeI,realEndDataI)
shading interp

%% Testing of stretched exponential fitting of these data
% addpath('/Users/tyler/Dropbox/Data/Frescoes with Kaori/Fragment1')
% 
% clear('xfit','ypred','coeffs','coeffserr','residuals')
% % y0_guess = 0;
% A_guess = 12e4;
% t2_guess = 0.1;
% expParam = 1;
% 
% 
% guesses = [A_guess;t2_guess;expParam]; %y0_guess;
% 
% CI = 90; %desired confidence interval in percent
% 
% for i = 1:length(depthList)
%     try
%     [xfit,ypred(:,i),coeffs(:,i),coeffserr(:,i),residuals(:,i)] = monodecay_t2_stretchedfit(echoTime,realEndData(:,i),guesses,CI);
%     catch
%         ypred(:,i) = zeros(length(xfit),1);
%         coeffs(:,i) = [0;0;0];
%         coeffserr(:,i) = [0;0;0];
%         residuals(:,i) = zeros(length(echoTime),1);
%     end
% end
% %%
% % xlims  =[-200 12000];
% xlims = [1000 5000];
% posnOffset = 0;
% position = depthList;
% 
% figure(4)
% subplot(3,1,1)
% hold on
% plot(position-posnOffset,coeffs(1,:),'-b')
% plot(position-posnOffset,coeffs(1,:)-coeffserr(1,:),':b')
% plot(position-posnOffset,coeffs(1,:)+coeffserr(1,:),':b')
% % plot(position-posnOffset,real(cplxPhased(1,:)),'-g')
% % ylim([0-50000 50000])
% xlim(xlims)
% ylabel('signal intensity (arb)')
% subplot(3,1,2)
% hold on
% plot(position-posnOffset,coeffs(2,:),'-r')
% plot(position-posnOffset,coeffs(2,:)-coeffserr(2,:),':r')
% plot(position-posnOffset,coeffs(2,:)+coeffserr(2,:),':r')
% ylim([0 0.1])
% xlim(xlims)
% ylabel('T_2 time/ms')
% subplot(3,1,3)
% hold on
% plot(position-posnOffset,coeffs(3,:),'-k')
% plot(position-posnOffset,coeffs(3,:)-coeffserr(3,:),':k')
% plot(position-posnOffset,coeffs(3,:)+coeffserr(3,:),':k')
% ylim([0 1])
% xlim(xlims)
% ylabel('stretched exponential parameter alpha (arb)')
% xlabel('position/um') 
%% Do ILT
alpha = 1e7;
lowLim = min(echoTimeI);
hiLim = max(echoTimeI)*10;
nrILTSteps = length(echoTimeI);
clear('spectrum','tau','chisq')
spectrum = zeros(nrILTSteps,length(depthList));
tau = spectrum;
chisq = zeros(length(depthList),1);

tic
parfor i = 1:length(depthList)
%     tic
    [spectrum(:,i),tau(:,i),chisq(i),~     ] = upnnlsmooth1D(realEndDataI(:,i),echoTimeI,  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
    disp(strcat('No. ',num2str(i)));
%     toc
end
totalTime = toc;
disp(strcat('Time per loop: ',num2str(totalTime/length(depthList)),' sec.'));

save('results02.mat')
% ILTArea = sum(spectrum,1);
 beep;

%% Plot
figure(4)
subplot(4,4,[1 2 3 5 6 7 9 10 11])
surf(depthList,log10(tau),spectrum);
view([0 90])
colormap(hot)
caxis([0 1e4])
 shading interp
 set(gca,'YDir','Reverse')
 xlabel('position (um)')
 ylabel('log(T_2/ms)')
 zlabel('intensity')
 ylim([log10(lowLim) log10(hiLim)])
 xlim([min(depthList) max(depthList)])
subplot(4,4,[13 14 15])
plot(depthList,sum(realEndData));
 set(gca,'YDir','Reverse')
    xlabel('distance (um) (toward MOUSE <-- --> away from MOUSE)')
    xlim([min(depthList) max(depthList)])
subplot(4,4,[4 8 12])
plot(sum(spectrum'),log10(tau))
 ylim([log10(lowLim) log10(hiLim)])
 line([0 5000],[log10(echoTime(1)) log10(echoTime(1))])
 line([0 5000],[log10(echoTime(length(echoTime))) log10(echoTime(length(echoTime)))],'LineStyle','-.')
  set(gca,'YDir','Reverse')
%    ylabel('log(T_2/ms)')
