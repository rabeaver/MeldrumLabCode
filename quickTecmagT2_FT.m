clear
clc
close all

%%
filename = 'WeirdStuff.tnt';
filedir = 'C:\Users\jnking01\Desktop\';
fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);

zf = 2;

tEcho = 200; %us
nEchoes = 1024;
nPts = 56;
nPtsBlank = 0;
tD = 2e-6; %dwell time, s

G = 6.59;
gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1

%% SNR calc


[~,Spoint] = max(real(spec2));
Spoint = Spoint + 3.5*nPts;
S = (real(spec(Spoint-nPts/2:Spoint+nPts/2)));
N = (imag(spec(Spoint-nPts/2:Spoint+nPts/2)));
% N = (real(specN(Spoint-nPts/2:Spoint+nPts/2)))';

SNR = snr(S,N)

figure(1)
hold on
plot(S)
plot(N)

%%


echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

data = reshape(spec,nPts,nEchoes);
data = data(1:(nPts-nPtsBlank),:);
dataInt = sum(data,1);
dataIntRe = real(dataInt)./max(real(dataInt));
dataIntIm = imag(dataInt)./max(real(dataInt));

Fs = 1/tD;                           % Sampling frequency 
L = (nPts)*(2^zf);          % Length of signal
NFFT = 2^nextpow2(L);               % Next power of 2 from length of y

t = (-(L-1)/2:L/2)*tD;               % Time vector
f = linspace(-Fs/2,Fs/2,NFFT);      % Hz
z = f/(gamma*G);                    % um, 280.47 Hz/um (for PM25)

dat = padarray(data, size(data(:,1),1)/2*((2^zf)-1),0); % Pad with 0's

profiles = flipud(fftshift(fft(dat,NFFT)/L, 1)); % Performs FFT algorithm

close
plot(z,abs(profiles(:,1)));

%%
guess = [1 15e-3];% 0.6 6e-03];
beta = nlinfit(echoVector,dataIntRe, @t2monofit_simple, guess);
ypred = t2monofit_simple(beta,echoVector);

figure(2)
hold on
plot(echoVector,dataIntRe);
plot(echoVector,dataIntIm);
plot(echoVector,ypred,'-r');
xlabel('time[s]');
ylabel('intensity');
set(gca, 'Fontsize',18,'linewidth',2);