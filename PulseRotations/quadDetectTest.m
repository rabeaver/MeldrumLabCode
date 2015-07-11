close all
clc
clear

v1 = -25; %Hz
% v2 = -30; %Hz
T  = 1;  %generic relaxation time T (s)
phi = 0; %phase

dw = 10e-3; %s
N = 2^10;    %sample points

w1 = v1*2*pi; %rad s-1
% w2 = v2*2*pi; %rad s-1
sampleF = 1/dw;
tAcq = dw*N;
SW = 0.5/dw;
dres = 1/tAcq;

t = (0:1:N-1)*dw;
f = (-0.5*N:1:0.5*N-1)/(N*dw);

FID = cos(t.*w1+phi).*exp(-t/T); %+cos(t.*w2+phi).*exp(-t/T);
Sminus = fft(FID);
% Sminus = fft(-FID);

v1 = 25; %Hz
% v2 = 30; %Hz
% T  = 1;  %generic relaxation time T (s)
% phi = 0; %phase
% 
% dw = 10e-3; %s
% N = 2^9;    %sample points
% 
w1 = v1*2*pi; %rad s-1
% w2 = v2*2*pi; %rad s-1
% sampleF = 1/dw;
% tAcq = dw*N;
% SW = 0.5/dw;
% dres = 1/tAcq;
% 
% t = (0:1:N-1)*dw;
% f = (-0.5*N:1:0.5*N-1)/(N*dw);

FID = cos(t.*w1+phi).*exp(-t/T); %+cos(t.*w2+phi).*exp(-t/T);
Splus = fft(FID);
% Sminus = fft(-FID);

%%
figure(1)
hold on
plot(t,FID)
title('Artificial FID, only real signal, 250 Hz')
xlabel('time [s]') 

figure(2)
hold on
plot(f,real(Splus),'-k')
plot(f,imag(Splus),'-r')
title('Artificial spectrum, 250 Hz')
xlabel('frequency [Hz]')
legend('Real','Imag')

%%
FIDi = sin(t.*w1+phi).*exp(-t/T); %+sin(t.*w2+phi).*exp(-t/T);

S_Re = Splus;
S_Im = fft(FIDi);

figure(3)
hold on
plot(f,imag(S_Im),'-r')
plot(f,real(S_Re),'-k')
legend('Im(imaginary)','Re(real)')
title('Artificial spectrum from complex signal (cos + sin), 250 Hz')
xlabel('frequency [Hz]')
%%

figure(4)
hold on
plot(f,real(S_Re)+imag(S_Im),'-k')
% plot(f,real(S_Re)-imag(S_Im),'-r')
title('Result of quadrature detection, 250 Hz')
legend('sum of im(Imag) + re(Real)')
xlabel('frequency [Hz]')