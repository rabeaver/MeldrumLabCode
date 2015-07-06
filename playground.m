clear
clc
close all

%%
T = 0.512;  %time (ms)
L = 512;    %nr cplx pts

dt = T/L;   %dwell time (ms)
Fs = 1/dt;  %sampling freq (kHz)
t = dt:dt:T; %time vector

x_axis = (t - T/2);
real_data = sin(x_axis)./x_axis;
imag_data = (sin(x_axis)./x_axis).*acoth(x_axis);
cplx_data = complex(real_data,real(imag_data));
cplx_data(L/2) = 1;


%%
NFFT = 2^nextpow2(L);
freq_space = linspace(-Fs/2,Fs/2,NFFT); %kHz
real_space = freq_space/852.51; %mm

lb = 0;
g_x = exp(-(t-(max(t)/2)).^2/(2*(1/lb)^2));
data = cplx_data.*g_x;
Y = fftshift(fft(data,NFFT)/L);
dY = diff(abs(Y))./diff(real_space);

figure(1)

subplot(1,3,1)
hold on
plot(t,real(cplx_data.*g_x),'-b')
plot(t,g_x*max(real(cplx_data)),'-.k')
xlabel('time (ms)')
xlim([0 T])

subplot(1,3,2)
plot(real_space,abs(Y));
xlabel('distance (mm) (toward MOUSE <-- --> away from MOUSE)')
% xlim([-0.1 0.1])

subplot(1,3,3)
plot(real_space(2:end),dY)
% xlim([-0.1 0.1])

resolution_um = 1000*range(real_space)/L


%%

x = -200:0.1:200;
real = sin(x)./x;
imag = sin(x-pi/2)./x;
figure
hold on
plot(x,real,'-k')
plot(x,imag,'-b')
%%
N = 256;
spatial_x = linspace(-10,10,N);
spatial_y = ones(1,length(spatial_x));
spatial_y(round(length(spatial_x)/2):end) = 0;

freq_y = fftshift(fft(spatial_y,N)/N);

figure
hold on
plot(spatial_x,real(freq_y),'-k')
plot(spatial_x,imag(freq_y),'-b')


testy = (sin(spatial_x)./spatial_x).*acoth(spatial_x);