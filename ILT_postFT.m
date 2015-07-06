
clear
clc
close all
% try
%     matlabpool(4)
% catch
% end

%%
% for i = 1:1:19;
%     dirnum = num2str(i);
%     dirhead = 'C:/Users/TMeldrum/Dropbox/Data from MobileMOUSE/cream/dblphantom/';
%     dir = strcat(dirhead,dirnum,'/');
%     cd(dir)

i = 1; %echo number for plotting
reduction = 1; %reduce sampling rate by this factor to decrease bandwidth 
lb = 25; %line broadening
zerofilling = 2; %factor for ZF (0 = no ZF, 1 = double the length, 2 = 4x the length, etc...)
sample_window = 50; %window for saving samples in um

real_data = load('Re_echoData.dat');
imag_data = load('Im_echoData.dat');

cplx_data = (padarray(complex(real_data,imag_data)',(size(real_data,2)/2)*((2^zerofilling)-1),0))';

parfilestem = 'acqu';

params.acqTime = readpar_Kea(strcat(parfilestem,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(parfilestem,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(parfilestem,'.par'),'nrScans');
params.rxPhase = readpar_Kea(strcat(parfilestem,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(parfilestem,'.par'),'rxGain');
params.nrPts = readpar_Kea(strcat(parfilestem,'.par'),'nrPnts');
params.repTime = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.repTime = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.b1Freq = readpar_Kea(strcat(parfilestem,'.par'),'b1Freq');
params.nrEchoes = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoes');
params.echoTime = readpar_Kea(strcat(parfilestem,'.par'),'echoTime');


T = params.acqTime*(2^zerofilling);     %total sample time (ms)
L = params.nrPts*(2^zerofilling)/reduction;       %number of sampled points
dt = T/L;               %dwell time (ms)
Fs = 1/dt;              %sampling freq (kHz)
t = dt:dt:T;            %time vector

echoVector = params.echoTime:params.echoTime:params.echoTime*params.nrEchoes;

% figure
% hold on
% plot(t,sum(real_data),'-k')
% plot(t,sum(imag_data),'-r')
% xlabel('time (milliseconds)')

%% Do FT
NFFT = 2^nextpow2(L);
freq_space = linspace(-Fs/2,Fs/2,NFFT); %kHz
real_space = freq_space/852.51; %mm

g_x = exp(-(t-(max(t)/2)).^2/(2*(1/lb)^2));
Y = zeros(params.nrEchoes,L);
dY = zeros(params.nrEchoes,L-1);

for j = 1:params.nrEchoes;
    data =  sum(cplx_data(j,:),1);
    data_dec = decimate(data,reduction).*g_x;
    Y(j,:) = fftshift(fft(data_dec,NFFT)/L);
    dY(j,:) = diff(abs(Y(j,:)))./diff(real_space);
end

%% Plotting

figure(1)
title('summed echo data')
subplot(1,3,1)
    hold on
    plot(t,imag(decimate(sum(cplx_data),reduction).*g_x),'-r')
    plot(t,real(decimate(sum(cplx_data),reduction).*g_x),'-b')
    plot(t,g_x*max(sum(real(cplx_data))),'-.k')
    xlabel('time (ms)')
    xlim([0 T])
subplot(1,3,2)
    plot(real_space,abs(sum(Y)));
    xlabel('distance (mm) (toward MOUSE <-- --> away from MOUSE)')
    xlim([-0.4 0.4])
subplot(1,3,3)
    plot(real_space(2:end),sum(dY))
    xlim([-0.4 0.4])
    
%%
figure(2)
title('individual echo data')
subplot(1,3,1)
    hold on
    plot(t,imag(decimate(cplx_data(i,:),reduction).*g_x),'-r')
    plot(t,real(decimate(cplx_data(i,:),reduction).*g_x),'-b')
    plot(t,g_x*max(real(cplx_data(i,:))),'-.k')
    xlabel('time (ms)')
    xlim([0 T])
subplot(1,3,2)
    plot(real_space,abs(Y(i,:)));
    xlabel('distance (mm) (toward MOUSE <-- --> away from MOUSE)')
    xlim([-0.4 0.4])
subplot(1,3,3)
    plot(real_space(2:end),dY(i,:))
    xlim([-0.4 0.4])
    
%%
figure(3)
surf(real_space,echoVector,abs(Y));
view([0 90])
xlabel('real space (mm) (toward MOUSE <-- --> away from MOUSE)')
xlim([-0.3 0.3])
ylabel('time (\mus)')
zlabel('signal intensity (arb)')

%% Do ILT
alpha = 1e6;
lowLim = min(echoVector)/1000;
hiLim = max(echoVector)/10;
nrILTSteps = length(echoVector);
clear('spectrum','tau','chisq')
spectrum = zeros(nrILTSteps,length(real_space));
tau = spectrum;
chisq = zeros(length(real_space),1);

tic
parfor i = 1:length(real_space)
%     tic
    [spectrum(:,i),tau(:,i),chisq(i),~     ] = upnnlsmooth1D(abs(Y(:,i)),echoVector'/1000,  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
    disp(strcat('No. ',num2str(i)));
%     toc
end
totalTime = toc;
disp(strcat('Time per loop: ',num2str(totalTime/length(real_space)),' sec.'));

ILTArea = sum(spectrum,1);
% rawArea = sum(real(phasedData),1);
% save('processed_data.mat')

% figure
% semilogx(tau,spectrum)

%% Plot
figure(1)
subplot(4,4,[1 2 3 5 6 7 9 10 11])
surf(real_space,log10(tau),spectrum);
view([0 90])
colormap(hot)
 shading interp
 set(gca,'YDir','Reverse')
 xlabel('position (mm)')
 ylabel('log(T_2/ms)')
 zlabel('intensity')
 ylim([log10(lowLim) log10(hiLim)])
 xlim([-0.2 0.2])
subplot(4,4,[13 14 15])
plot(real_space,abs(sum(Y)));
 set(gca,'YDir','Reverse')
    xlabel('distance (mm) (toward MOUSE <-- --> away from MOUSE)')
 xlim([-0.2 0.2])
subplot(4,4,[4 8 12])
plot(sum(spectrum'),log10(tau))
 ylim([log10(lowLim) log10(hiLim)])
  set(gca,'YDir','Reverse')
%    ylabel('log(T_2/ms)')
