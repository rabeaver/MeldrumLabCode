clear
close all
clc

%% User parameters

% Input data location and filestem, number of profiles, depths per profile,
% and echoes per depth
% addpath('C:\Users\TMeldrum\Dropbox\Coding\Matlab Processing\Scripts\')

cd('C:\Users\tkmeldrum\Desktop\Bead Pack\1500uMGd_BEADPACK_50u\1\')
datafilestem = 'data2';
acqufilestem = 'acqu';
omit_points = 0;
auto_phasing = 0;

params.numberProfiles = readpar_Kea(strcat(acqufilestem,'.par'),'nrExp');
params.finalDepth = readpar_Kea(strcat(acqufilestem,'.par'),'finalDepth');
params.initDepth = readpar_Kea(strcat(acqufilestem,'.par'),'initDepth');
params.stepSize = readpar_Kea(strcat(acqufilestem,'.par'),'stepSize');
params.numberDepths = 1 + (params.finalDepth - params.initDepth)/params.stepSize;
params.numberEchoes = readpar_Kea(strcat(acqufilestem,'.par'),'nrEchoes');
params.rxPhase = readpar_Kea(strcat(acqufilestem,'.par'),'rxPhase');
params.bandwidth = readpar_Kea(strcat(acqufilestem,'.par'),'bandwidth');
params.pulseLength = readpar_Kea(strcat(acqufilestem,'.par'),'pulseLength');
params.b1Freq = readpar_Kea(strcat(acqufilestem,'.par'),'b1Freq');
params.numberScans = readpar_Kea(strcat(acqufilestem,'.par'),'nrScans');
params.repTime = readpar_Kea(strcat(acqufilestem,'.par'),'repTime');
params.echoTime = readpar_Kea(strcat(acqufilestem,'.par'),'echoTime');

% Load data
data = dlmread(strcat(datafilestem,'.csv'));
echotime = data(:,1);
if auto_phasing == 1;
    signal = complex(data(:,2),data(:,3));
    amplitude = real(autophase(signal,0.1));
else
    amplitude = data(:,2:end);
end
% %% Plot
% 
figure(1)
hold on
scatter(echotime,real(signal),'ok')
scatter(echotime,imag(signal),'or')
xlabel('echo time (ms)')
xlim([0 max(echotime)])

%% Fit to exp decay

y0_guess = 1;
A_guess = max(amplitude(omit_points+1:end,1));
t2_guess = 5;
A_guess2 = A_guess;
t2_guess2 = 2; %ms

bounds = [-10 10;
          0 5000;
          0 100;
          0 5000;
          0 10];

guesses = [y0_guess;A_guess;t2_guess]; %;A_guess2;t2_guess2];

CI = 90; %desired confidence interval in percent

% for i = 1:1:size(amplitude,2)
    [xfit,ypred,coeffs,residuals] = monodecay_t2fit(echotime(omit_points+1:end),amplitude(omit_points+1:end,1),guesses,CI);
% end

%     [xfit,ypred,coeffs,residuals] = bidecay_t2fit(echotime(omit_points+1:end),amplitude(omit_points+1:end,1),guesses,CI);

fid = fopen('t2vals.txt','a+');
fprintf(fid,'%i %f %f %f %f %f %f\n', params.echoTime,coeffs(1,1), coeffs(1,2),coeffs(2,1), coeffs(2,2),coeffs(3,1), coeffs(3,2));
fclose(fid);
 
%%
% textinfo(5) = {'\tau_{t2}'};
figure(2)
subplot(2,1,1)
hold all
% scatter(echotime(omit_points+1:end),amplitude(omit_points+1:end,1),'or')
scatter(echotime,amplitude(:,1),'or')
plot(xfit,ypred)
xlabel('ms')
ylabel('signal amplitude')
subplot(2,1,2)
hold on
plot(echotime(omit_points+1:end),residuals,'-b')
plot(echotime(omit_points+1:end),zeros(length(residuals)),'-k')
% text(0.4*max(xfit),0.8*max(ypred),textinfo)
