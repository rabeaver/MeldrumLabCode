
% FOR USE WITH the Museum Ludwig data collected on original paintings, Jan
% 2012

clear
clc
close all
try
    matlabpool
catch
end

%%
% for i = 1:1:19;
%     dirnum = num2str(i);
%     dirhead = 'C:/Users/TMeldrum/Dropbox/Data from MobileMOUSE/cream/dblphantom/';
%     dir = strcat(dirhead,dirnum,'/');
%     cd(dir)

omitpoints = 1;
i = 1; %echo number for plotting
% reduction = 1; %reduce sampling rate by this factor to decrease bandwidth 
% lb = 25; %line broadening
% zerofilling = 2; %factor for ZF (0 = no ZF, 1 = double the length, 2 = 4x the length, etc...)
% sample_window = 50; %window for saving samples in um


filename = 'pechstein_haus_pos3_profile2';
real_data = load(strcat(filename,'-decaysRe.dat'));
imag_data = load(strcat(filename,'-decaysIm.dat'));

contrast_data = load(strcat(filename,'.dat'));
position = contrast_data(:,1);

cplx_data = complex(real_data(omitpoints+1:end,2:end),imag_data(omitpoints+1:end,2:end));

params.acqTime = readpar_Kea(strcat(filename,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(filename,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(filename,'.par'),'nrScans');
params.rxPhase = readpar_Kea(strcat(filename,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(filename,'.par'),'rxGain');
params.nrPts = readpar_Kea(strcat(filename,'.par'),'nrPnts');
params.repTime = readpar_Kea(strcat(filename,'.par'),'repTime');
params.repTime = readpar_Kea(strcat(filename,'.par'),'repTime');
params.b1Freq = readpar_Kea(strcat(filename,'.par'),'b1Freq');
params.nrEchoes = readpar_Kea(strcat(filename,'.par'),'nrEchoes');
params.echoTime = readpar_Kea(strcat(filename,'.par'),'echoTime');


% T = params.acqTime*(2^zerofilling);     %total sample time (ms)
% L = params.nrPts*(2^zerofilling)/reduction;       %number of sampled points
% dt = T/L;               %dwell time (ms)
% Fs = 1/dt;              %sampling freq (kHz)
% t = dt:dt:T;            %time vector

echoVector = params.echoTime*(omitpoints+1):params.echoTime:params.echoTime*params.nrEchoes;
gaussianWeight = exp(-echoVector/(max(echoVector)));
phased_data = autophase(cplx_data,1);
for i = 1:length(position)
    weighted_data(:,i) = phased_data(:,i)'.*gaussianWeight;
end
% zfArray = padarray(weighted_data,size(weighted_data,1),'post');
% zfEchoVector = params.echoTime:params.echoTime:2*params.echoTime*params.nrEchoes;

% figure
% hold on
% plot(t,sum(real_data),'-k')
% plot(t,sum(imag_data),'-r')
% xlabel('time (milliseconds)')

%% Plotting

figure(1)
plot(position,real(cplx_data(1,:)))
set(gca,'XDir','Reverse')

figure(2)
plot(position,real(phased_data(1,:)))
set(gca,'XDir','Reverse')

figure(3)
surf(position,echoVector,real(cplx_data))
shading interp


%% Do ILT
alpha = 2e7;
lowLim = min(echoVector)/1000;
hiLim = max(echoVector)/100;
nrILTSteps = length(echoVector);
clear('spectrum','tau','chisq')
spectrum = zeros(nrILTSteps,length(position));
tau = spectrum;
chisq = zeros(length(position),1);

tic
parfor i = 1:length(position)
%     tic
    [spectrum(:,i),tau(:,i),chisq(i),~     ] = upnnlsmooth1D(real(cplx_data(:,i)),echoVector'/1000,  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
    disp(strcat('No. ',num2str(i)));
%     toc
end
totalTime = toc;
disp(strcat('Time per loop: ',num2str(totalTime/length(position)),' sec.'));

ILTArea = sum(spectrum,1);
beep;

% rawArea = sum(real(phasedData),1);
% save('processed_data.mat')

% figure
% semilogx(tau,spectrum)

%% Plot
figure(4)
subplot(4,4,[1 2 3 5 6 7 9 10 11])
surf(position,log10(tau),spectrum);
view([0 90])
colormap(hot)
caxis([0 30])
 shading flat
 set(gca,'YDir','Reverse')
 xlabel('position (mm)')
 ylabel('log(T_2/ms)')
 zlabel('intensity')
 ylim([log10(lowLim) log10(hiLim)])
%  xlim([1800 2800])
subplot(4,4,[13 14 15])
plot(position,sum(real(cplx_data)));
 set(gca,'YDir','Reverse')
    xlabel('distance (mm) (toward MOUSE <-- --> away from MOUSE)')
%      xlim([1800 2800])
subplot(4,4,[4 8 12])
plot(sum(spectrum'),log10(tau))
 ylim([log10(lowLim) log10(hiLim)])
 line([0 5000],[log10(echoVector(1)/1000) log10(echoVector(1)/1000)])
 line([0 5000],[log10(echoVector(length(echoVector))/1000) log10(echoVector(length(echoVector))/1000)],'LineStyle','-.')
  set(gca,'YDir','Reverse')
%    ylabel('log(T_2/ms)')

%%
figure(5)
subplot(4,4,[1 2 3 5 6 7 9 10 11])
hold on
surf(position4,log10(tau4),spectrum4);
view([0 90])
colormap(cool)
 shading interp
 set(gca,'YDir','Reverse')
 xlabel('position (mm)')
 ylabel('log(T_2/ms)')
 zlabel('intensity')
 ylim([log10(lowLim) log10(hiLim)])
 xlim([1800 2800])
subplot(4,4,[13 14 15])
plot(position4,sum(real(phased_data4)));
 set(gca,'YDir','Reverse')
    xlabel('distance (mm) (toward MOUSE <-- --> away from MOUSE)')
    xlim([1800 2800])
subplot(4,4,[4 8 12])
plot(sum(spectrum4'),log10(tau4))
 ylim([log10(lowLim) log10(hiLim)])
  set(gca,'YDir','Reverse')
%    ylabel('log(T_2/ms)')