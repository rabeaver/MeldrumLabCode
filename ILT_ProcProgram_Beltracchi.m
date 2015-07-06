
% FOR USE WITH the Beltracchi data, collected Feb 2012 (LKA Berlin)
clear
clc
close all
try
    matlabpool
catch
end

%%
omitpoints = 1;

positionList = 'depthInfo.dat';
parname = 'acqu';

positionInfo = load(positionList);
allPositions = positionInfo(:,2);
positionCatch = positionInfo(:,4);
position = allPositions.*positionCatch;
position = position(position~=0);

params.acqTime = readpar_Kea(strcat(parname,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(parname,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(parname,'.par'),'nrScansT2');
params.rxPhase = readpar_Kea(strcat(parname,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(parname,'.par'),'rxGain');
params.nrPts = readpar_Kea(strcat(parname,'.par'),'nrPnts');
params.repTime = readpar_Kea(strcat(parname,'.par'),'repTimeT2');
params.b1Freq = readpar_Kea(strcat(parname,'.par'),'b1Freq');
params.nrEchoes = readpar_Kea(strcat(parname,'.par'),'nrEchoesT2');
params.echoTime = readpar_Kea(strcat(parname,'.par'),'echoTime');

echoVector = params.echoTime*(omitpoints+1):params.echoTime:params.echoTime*params.nrEchoes;

real_data = zeros(length(echoVector),length(position));
imag_data = real_data;

for i = 1:length(position)
    temp = load(strcat('T2data-',num2str(position(i)),'.dat'));
    real_data(:,i) = temp(omitpoints+1:end,2);
    imag_data(:,i) = temp(omitpoints+1:end,3);
end

cplx_data = complex(real_data,imag_data);


%% Plotting

figure(1)
plot(position,real(cplx_data(1,:)))
set(gca,'XDir','Reverse')

figure(2)
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
    [spectrum(:,i),tau(:,i),chisq(i),~     ] = upnnlsmooth1D(real(cplx_data(:,i)),echoVector'/1000,  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
    disp(strcat('No. ',num2str(i)));
end
totalTime = toc;
disp(strcat('Time per loop: ',num2str(totalTime/length(position)),' sec.'));

ILTArea = sum(spectrum,1);
beep;


%% Plot
figure(4)
subplot(4,4,[1 2 3 5 6 7 9 10 11])
surf(position,log10(tau),spectrum);
view([0 90])
colormap(hot)
caxis([0 1e4])
 shading flat
 set(gca,'YDir','Reverse')
 xlabel('position (mm)')
 ylabel('log(T_2/ms)')
 zlabel('intensity')
%  ylim([log10(lowLim) log10(hiLim)])
%  xlim([1800 2800])
   zlim([0 1e4])
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
% figure(5)
% subplot(4,4,[1 2 3 5 6 7 9 10 11])
% hold on
% surf(position4,log10(tau4),spectrum4);
% view([0 90])
% colormap(cool)
%  shading interp
%  set(gca,'YDir','Reverse')
%  xlabel('position (mm)')
%  ylabel('log(T_2/ms)')
%  zlabel('intensity')
%  ylim([log10(lowLim) log10(hiLim)])
%  xlim([1800 2800])
% subplot(4,4,[13 14 15])
% plot(position4,sum(real(phased_data4)));
%  set(gca,'YDir','Reverse')
%     xlabel('distance (mm) (toward MOUSE <-- --> away from MOUSE)')
%     xlim([1800 2800])
% subplot(4,4,[4 8 12])
% plot(sum(spectrum4'),log10(tau4))
%  ylim([log10(lowLim) log10(hiLim)])
%   set(gca,'YDir','Reverse')
% %    ylabel('log(T_2/ms)')