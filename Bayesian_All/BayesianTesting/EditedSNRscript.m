clear
clc
close all
%% Load Tecmag Data
[params,spec1,spec2] = readTecmag('/Users/tyler/Desktop/NAU/water_noiseCPMGnoise_n256_12June2014.tnt');
nrPnts = 69;
nrEchoes = 320;
nrCalcEchoes = 320;
nrBlankPoints = 5;
realPts = nrPnts - nrBlankPoints;
echoVector = linspace(1,100,nrEchoes);
nrNoisePoints = 1e6;

spec3 = reshape(spec1',nrPnts,nrEchoes); %move from one vector for all echoes to each echo separately
spec3 = spec3(1:end-nrBlankPoints,1:nrCalcEchoes,:);
t1 = (1:realPts)'*params.dw(1); %time axis for eches
spec2 = sum(spec3,3);
spec3 = spec2(1:end,129:192);
% specN = [spec2(1:end,1:128) spec2(1:end,193:320)];
% specN = reshape(specN,1,64*256);
%params for echo shape determination
v(1) = 7e-3;     %psi
v(2) = 1e-5;    %sigma
v(3) = 2e6;     %amp
v(4) = 42e-6;   %echo center (s)

% the echo shape function is now in 'leastsquares.m'
% hfxn = @(v,t)v(3)*exp(-2*pi*1i*v(1).*(t-v(4))).*exp(-(t-v(4)).^2./(2*v(2)^2));
opts = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','Display','off'); %fitting options

%% Determine Echo Shape
opts = optimoptions(@lsqnonlin,opts); %fitting options
phat = zeros(64,4); %preallocate
for i = 1:64
    phat(i,:) = lsqnonlin('leastsquares',v,[], [], opts, t1, spec3(:,i));
end
phat(:,2) = abs(phat(:,2));
meanphat = mean(phat,1); %avg value for p-hat (all echo shape parameters)
meanpsi = meanphat(1);
stdpsi = std(phat(:,1));
meansigma = meanphat(2);
stdsigma = std(phat(:,2));
meanamp = meanphat(3);
stdamp = std(phat(:,3));
meancenter = meanphat(4);
stdcenter = std(phat(:,4));

%% Sample figure for checking
figure(1) %plot a few echoes and their calculated shape for comparison
for i=1:16
subplot(4,4,i)
hold on
plot(t1,real(CplxGaussian(phat(i,:),t1)),'-k')
plot(t1,imag(CplxGaussian(phat(i,:),t1)),'-r')
plot(t1,real(spec3(:,i,1)),'-k')
plot(t1,imag(spec3(:,i,1)),'-r')
plot(t1,abs(spec3(:,i,1)),'-b')
plot(t1,abs(CplxGaussian(phat(i,:),t1)),'-b')
end

%% Make echo shape vector
stdEchoParams = [meanpsi, meansigma, 1, nrPnts*params.dw(1)/2];
% stdEchoParams = [meanpsi, meansigma, 1, meancenter];
h = [real(CplxGaussian(stdEchoParams,t1));imag(CplxGaussian(stdEchoParams,t1))];

%% Make fake noise
%generate fake noise using points 1-20 of echo 5 of set 1 (spec3(1:20,5,1))
%until can get the real noise generator working. I match the mean and
%standard deviation so the noise is properly scaled. 
% xx=10;
% close all
% 
% stdNoise = spec3(1:20,5);
% noiseParams = [ mean(real(stdNoise)) std(real(stdNoise));
%                 mean(imag(stdNoise)) std(imag(stdNoise))];

%two sets of fake noise for before and after signal acquisition
% fakeNoise = complex(randn(nrNoisePoints,1)*noiseParams(1,2)+noiseParams(1,1),randn(nrNoisePoints,1)*noiseParams(2,2)+noiseParams(2,1));
[ap,noise,~] = readTecmag('/Users/tyler/Desktop/NAU/eraser_noise256.tnt');

% noise = [spec2(:,1:128) spec2(:,193:320)];
% noise = reshape(noise,256*64,1);

% spec = spec2(
% realNoise = spec3(:,22:32);
% fakeNoise = reshape(realNoise,59*11,1);
%plot to check the scaling of the noise
% figure(2)
% hold on
% plot(t1,real(spec3(:,5,1)),'-k')
% plot(t1,imag(spec3(:,5,1)),'-r')
% plot(t1,real(CplxGaussian(phat(5,:),t1)),'-k')
% plot(t1,imag(CplxGaussian(phat(5,:),t1)),'-r')
% plot(t1,real(fakeNoise),'-b')
% plot(t1,imag(fakeNoise),'-g')
% plot(t1,real(fakeNoise2),'-b')
% plot(t1,imag(fakeNoise2),'-g')

% Noise covariance
%this section makes me nervous--I don't really understand what's going on,
%but the dimensions agree and it's a cross(/auto)covariance

% [c11,lags11] = xcov(real(fakeNoise),real(fakeNoise'),nrPnts-nrBlankPoints);
[c11,lags11] = xcov(real(noise),real(noise),realPts,'coeff');
[c22,lags22] = xcov(imag(noise),imag(noise),realPts,'coeff');
[c12,lags12] = xcov(real(noise),imag(noise),realPts,'coeff');
[c21,lags21] = xcov(imag(noise),real(noise),realPts,'coeff');
% id = ones(realPts);
tildeRange = [ceil(2*realPts-1.5*realPts) floor(2*realPts-0.5*realPts)-1];
% tildeRange = [65 128];
Ctilde = [diag(c11(tildeRange(1):tildeRange(2))), diag(c12(tildeRange(1):tildeRange(2))); diag(c21(tildeRange(1):tildeRange(2))), diag(c22(tildeRange(1):tildeRange(2)))];
% Cnew = randn(118,118)*Ctilde;
Wtilde = inv(Ctilde);



%%
ahat_1 = inv(h'*Wtilde*h)*h'*Wtilde;
Y = [real(spec3);imag(spec3)];
ahat = ahat_1*Y;
%% plot(ahat)


data = h*ahat;
data = complex(data(1:nrPnts-nrBlankPoints,:),data(nrPnts-nrBlankPoints+1:end,:));
figure(2)
subplot(1,2,1)
surf(real(data))
shading flat
subplot(1,2,2)
surf(real(spec3(:,1:64)))
shading flat

%
intDataCalc = max(max(abs(data),1));
intDataMeas = sum(abs(spec3(:,1:64)),1);
figure(3)
hold on
plot(echoVector(1:64),intDataCalc,'-b')
plot(echoVector(1:64),intDataMeas,'-k')
legend('Calc','Measured')

figure(4)
hold on
plot(lags11,c11,'-k')
plot(lags22,c22,'-r')
plot(lags12,c12-0.5,'--k')
plot(lags21,c21-0.5,'--r')
legend('11','22','12','21')

figure(6)
hold on
plot(echoVector(1:64),intDataCalc./max(intDataCalc),'-b')
plot(echoVector(1:64),intDataMeas./max(intDataMeas),'-k')
legend('Calc','Measured')
% cftool(echoVector(1:64),intDataCalc)