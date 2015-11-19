% SNR Calculator

surf(data2d); shading flat
% Get peak signal and RMS in noise range

figure()
plot(data2d(79,:))

%Choose point
peak = 79; %Index of T1 point w/ highest amplitude
echo = 1; %Echo containing highest amp

%% Sum Noise
close all

figure()
plot(data2d(peak,:));

numNoise = 4; %Number of echoes determined to be noise

figure()
plot(data2d(peak,end-numNoise:end));

Noise = data2d(peak,end-numNoise:end);

RMSnoise = rms(Noise);


%% Calculate SNR

SNR = data2d(peak,echo)/RMSnoise;