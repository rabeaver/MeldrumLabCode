clear
clc
close all

%%
% for i = 1:1:19;
%     dirnum = num2str(i);
%     dirhead = 'C:/Users/TMeldrum/Dropbox/Data from MobileMOUSE/cream/dblphantom/';
%     dir = strcat(dirhead,dirnum,'/');
%     cd(dir)

i = 2; %echo number for plotting
reduction = 1; %reduce sampling rate by this factor to decrease bandwidth 
lb = 0; %line broadening
zerofilling = 0; %factor for ZF (0 = no ZF, 1 = double the length, 2 = 4x the length, etc...)
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
meshc(real_space,echoVector,abs(Y));
xlabel('real space (mm) (toward MOUSE <-- --> away from MOUSE)')
xlim([-0.3 0.3])
ylabel('time (\mus)')
zlabel('signal intensity (arb)')
%%
figure(4)
plot(real_space,abs(Y(i,:)))
xlabel('real space (mm) (toward MOUSE <-- --> away from MOUSE)')
xlim([-0.3 0.3])

% figure(5)
% plot(real_space(2:end)-(real_space(2)-real_space(1))/2,diff(abs(Y(i,:))));
% xlabel('real space (mm) (toward MOUSE <-- --> away from MOUSE)')
% xlim([-0.15 0.15])

resolution_um = 1000*reduction*range(real_space)/NFFT;
samplerange = round(sample_window/resolution_um/2);
profilerange = abs(Y(1,length(Y)/2 - samplerange:length(Y)/2 + samplerange-1))';
profilex = (1:2*samplerange)'*resolution_um + (i-1)*(resolution_um*2*samplerange);

%% Define limits for T2 fitting (stacking multiple T2 curves together to get one improved signal)
fitlims = [31 38];

figure(6)
hold on
plot(echoVector,abs(Y(:,fitlims(1):fitlims(2))));

fitdata = [echoVector',abs(Y(:,fitlims(1):fitlims(2)))];
fitposition = real_space(fitlims(1):fitlims(2))';

masterfit_time = repmat(echoVector',length(fitposition),1);
masterfit_data = reshape(fitdata(:,2:end),length(fitposition)*length(echoVector),1);
masterfit = [masterfit_time, masterfit_data];
masterfit = sortrows(masterfit);

%% Do T2 fitting
y0_guess = 1;
A_guess = max(fitdata(1,2:end));
t2_guess = 2000;
% A_guess2 = A_guess;
% t2_guess2 = 2; %ms

guesses = [y0_guess;A_guess;t2_guess]; %;A_guess2;t2_guess2];

CI = 90; %desired confidence interval in percent

fid = fopen('t2vals.txt','w');
fprintf(fid,'position, y0, y0_err, A, A_err, T2, T2_err\n');
fclose(fid);

% for i = 1:size(fitdata,2)-1
%     [xfit,ypred,coeffs,residuals] = monodecay_t2fit(fitdata(:,1),fitdata(:,i+1),guesses,CI);
%     fid = fopen('t2vals.txt','a+');
%         fprintf(fid,'%f, %f, %f, %f, %f, %f, %f\n', fitposition(i), coeffs(1,1), coeffs(1,2), coeffs(2,1), coeffs(2,2),coeffs(3,1), coeffs(3,2));
%     fclose(fid);
% end
%     [xfit,ypred,coeffs,residuals] = monodecay_t2fit(masterfit(:,1),masterfit(:,2),guesses,CI);
%     fid = fopen('t2vals.txt','a+');
%         fprintf(fid,'%s, %f, %f, %f, %f, %f, %f\n', 'All', coeffs(1,1), coeffs(1,2), coeffs(2,1), coeffs(2,2),coeffs(3,1), coeffs(3,2));
%     fclose(fid);
%     [xfit,ypred,coeffs,residuals] = bidecay_t2fit(echotime(omit_points+1:end),amplitude(omit_points+1:end,1),guesses,CI);


%
% f = Fs/2*linspace(0,1,NFFT/2+1);
% f = Fs/2*linspace(-1,1,NFFT);

















%%
fid = fopen(strcat(dirhead,'profilex.txt'),'a+');
fprintf(fid,'%4.4f \n', profilex);
fclose(fid);
fid = fopen(strcat(dirhead,'profiley.txt'),'a+');
fprintf(fid,'%4.4f \n', profilerange);
fclose(fid);

% if i ~= 19
%     clear
%     close all
% end

% end
%%
cd(dirhead)
profile = load('profile.txt');
space = resolution_um:resolution_um:length(profile)*resolution_um;
scatter(space,profile)
%%
xdata = load('profilex.txt');
ydata = load('profiley.txt');
figure
plot(xdata,ydata)
xlabel('position (\mum)')
ylabel('signal intensity (arb)')

%%
g_x = exp(-(t-(NFFT/2*dt)).^2/(2*5e-3^2));