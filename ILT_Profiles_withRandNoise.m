
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

i = 1; %echo number for plotting
% reduction = 1; %reduce sampling rate by this factor to decrease bandwidth 
% lb = 25; %line broadening
% zerofilling = 2; %factor for ZF (0 = no ZF, 1 = double the length, 2 = 4x the length, etc...)
% sample_window = 50; %window for saving samples in um


filename = 'sample4_test6';
nrNoises = 1;

real_data = load(strcat(filename,'-decaysRe.dat'));
imag_data = load(strcat(filename,'-decaysIm.dat'));

contrast_data = load(strcat(filename,'.dat'));
position = contrast_data(:,1);

cplx_data = complex(real_data(:,2:end),imag_data(:,2:end));

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


echoVector = params.echoTime:params.echoTime:params.echoTime*params.nrEchoes;

noiseData = zeros(size(cplx_data,1),size(cplx_data,2),nrNoises);

parfor i=1:nrNoises;
    noiseData(:,:,i) = cplx_data + 0.05*max(max(real(cplx_data)))*randn([size(cplx_data,1),size(cplx_data,2)]);
end

%% Plotting

figure(1)
plot(position,real(cplx_data(1,:)))
set(gca,'XDir','Reverse')

figure(2)
surf(position,echoVector,real(cplx_data))
shading interp


%% Do ILT
alpha = 1e7;
lowLim = min(echoVector)/10000;
hiLim = max(echoVector);
nrILTSteps = length(echoVector);
clear('spectrum','tau','chisq')
spectrum = zeros(nrILTSteps,length(position));
tau = spectrum;
chisq = zeros(length(position));

tic

    parfor i = 1:length(position)
        %     tic
        [spectrum(:,i),tau(:,i),chisq(i),~     ] = upnnlsmooth1D(real(cplx_data(:,i)),echoVector'/1000,  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
        disp(strcat('No. ',num2str(i)));
        %     toc
    end

totalTime = toc;
disp(strcat('Time per loop: ',num2str(totalTime/length(position)),' sec.'));

% ILTArea = sum(spectrum,1);
beep;
%%
% close all
% 
% for i = 1:nrNoises;
%     figure(i)
%     surf(spectrum(:,:,i))
%     view([0 90])
% end
% 
spectrumAvg = mean(spectrum,3);
spectrumStd = std(spectrum,0,3);
spectrumVar = var(spectrum,0,3);
% 
% figure(1)
% surf(spectrumAvg)
% view([0 90])

% rawArea = sum(real(phasedData),1);
% save('processed_data.mat')

% figure
% semilogx(tau,spectrum)

%% Plot
figure(4)
subplot(4,4,[1 2 3 5 6 7 9 10 11])
contourf(position,log10(tau(:,1,1)),spectrum,100);
view([0 90])
colormap(hot)
caxis([0 50])
 shading flat
 set(gca,'YDir','Reverse')
 xlabel('position (mm)')
 ylabel('log(T_2/ms)')
 zlabel('intensity')
 ylim([log10(echoVector(1)/1000) log10(hiLim)])
 zlim([0 50])
%  xlim([1800 2800])
subplot(4,4,[13 14 15])
plot(position,sum(real(cplx_data)));
 set(gca,'YDir','Reverse')
    xlabel('distance (mm) (toward MOUSE <-- --> away from MOUSE)')
%      xlim([1800 2800])
subplot(4,4,[4 8 12])
plot(sum(spectrum'),log10(tau(:,:,1)))
 ylim([log10(echoVector(1)/1000) log10(hiLim)])
%  line([0 5000],[log10(echoVector(1)/1000) log10(echoVector(1)/1000)])
 line([0 5000],[log10(echoVector(length(echoVector))/1000) log10(echoVector(length(echoVector))/1000)],'LineStyle','-.')
  set(gca,'YDir','Reverse')
%    ylabel('log(T_2/ms)')

%% Plot
figure(5)
subplot(4,4,[1 2 3 5 6 7 9 10 11])
contourf(position,log10(tau(:,1,1)),spectrumAvg,100);
view([0 90])
colormap(hot)
caxis([0 50])
 shading flat
 set(gca,'YDir','Reverse')
 xlabel('position (mm)')
 ylabel('log(T_2/ms)')
 zlabel('intensity')
 ylim([log10(echoVector(1)/1000) log10(hiLim)])
  line([0 5000],[log10(echoVector(length(echoVector))/1000) log10(echoVector(length(echoVector))/1000)],'LineStyle','-.')
% ylim([0 1])
 zlim([0 50])
%  xlim([1800 2800])
subplot(4,4,[13 14 15])
plot(position,sum(real(cplx_data)));
 set(gca,'YDir','Reverse')
    xlabel('distance (mm) (toward MOUSE <-- --> away from MOUSE)')
%      xlim([1800 2800])
subplot(4,4,[4 8 12])
plot(sum(spectrumAvg'),log(tau(:,1,1)))
%  ylim([log10(lowLim) log10(hiLim)])
ylim([log10(echoVector(1)/1000) log10(hiLim)])
% ylim([0 1])
%  line([0 5000],[log10(echoVector(1)/1000) log10(echoVector(1)/1000)])
 line([0 5000],[log10(echoVector(length(echoVector))/1000) log10(echoVector(length(echoVector))/1000)],'LineStyle','-.')
  set(gca,'YDir','Reverse')
%    ylabel('log(T_2/ms)')
