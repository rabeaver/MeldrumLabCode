clear
clc
close all

%% User parameters
filename = 'EtGlyLarge_STEr_14Jan2016.tnt';
filedir = 'C:\CommonData\EthyleneGlycol\';
fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tEcho = 280; %us
nEchoes = 128;
nPts = 50;
nPtsBlank = 0;
DELTA = 1e-3; %s
deltamin = 100e-6; %s
deltamax = 2600e-6; %s
refocused_delta = 1; % if there are 180 pulses in the delta periods, set to 1 to adjust the delta time
GM = 280; %MHz m-1


%% Non-user parameters
gamma = 2.675222005e8; %s-1 T-1
gammaM = 42.57748; %MHz T-1
G = GM/gammaM; %T m-1

n2DPts = ap.td(2);
echoVector = (tEcho:tEcho:nEchoes*tEcho)*1e-6;

if refocused_delta == 1
    delta = linspace(deltamin,deltamax,n2DPts)/2;
else
    delta = linspace(deltamin,deltamax,n2DPts);
end

data = reshape(spec2',nPts,nEchoes,n2DPts);
data = data(1:(nPts-nPtsBlank),:,:);
dataInt = sum(sum(data,1),2);
dataInt = reshape(dataInt,1,n2DPts);
% dataInt = dataInt./-dataInt(1);
dataIntRe = real(dataInt);
dataIntIm = imag(dataInt);

s1 = log10(dataIntRe./max(dataIntRe));
s_x = -gamma^2*G^2.*delta.^2.*(DELTA+delta./3)*1e-9;

%% Fit and give D
[p,S] = polyfit(s_x,s1,1);
D = p(1)*1e-9 %output value of D in m2/s

figure(1)
hold on
scatter(s_x,s1)
plot(s_x,polyval(p,s_x));
set(gca,'defaulttextinterpreter','latex')
xlabel('$-\gamma^{2}G^{2}\delta^{2}(\Delta+\frac{\delta}{3})\times 10^{-9}$')
ylabel('$log \frac{S}{S_0}$')

excelOut = [s1',s_x'*1e9];