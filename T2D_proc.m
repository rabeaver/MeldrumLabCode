clear
clc
close all

%%
datadir = 'C:\CommonData\EthyleneGLycol\';
datafile = 'EtGlyLarge_STE_refoc180_15Jan2016_1';


nPts = 76;                          % # of acqu points
omitPts = 4;                        % the number of points that are zeros from the spectrometer
nEchoes = 128;                      % Echoes
omitEchoes = 0;                     % numner of echoes to remove from data
tD = 8e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 700;                           % us
deltaMin = 400e-6;                  % s
deltaMax = 1600e-6;                 % s
DELTA = 10e-3;                      % s
noisePoints = 8;                   % number of points for measuring noise
noiseNumber = 1;                    % scan number to use for determining SNR
G = 6.59;                           % T m-1, B0 field gradient



gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1


[ap , spec, spec2] = readTecmag4d(strcat(datadir,datafile,'.tnt'));
T2Ddat = reshape(spec2, ap.td(2), nPts, nEchoes);
T2Ddat = T2Ddat(:,1:nPts-omitPts,omitEchoes+1:end);

deltaVec = linspace(deltaMin,deltaMax,ap.td(2));
xD = -gammaRad^2*G^2.*deltaVec.^2.*(DELTA+deltaVec/3)*1e-9;

echoVec = tE*(omitEchoes+1):tE:(nEchoes*tE);

data = sum(real(T2Ddat),2);
data = reshape(data,ap.td(2),(nEchoes-omitEchoes));


for i = 1:ap.td(2)
    dataY(i,:) = log10(abs(data(i,:)./data(i,1)));
end
%% SNR calc 
n1 = T2Ddat(noiseNumber,1:noisePoints,:);
n2 = T2Ddat(noiseNumber,nPts-noisePoints-omitPts:end,:);
n = cat(2,n1,n2);
n = reshape(n,1,(2*noisePoints+1)*(nEchoes-omitEchoes));
s = T2Ddat(noiseNumber,:,:);
s = reshape(s,1,(nPts-omitPts)*(nEchoes-omitEchoes));

figure
hold on
plot(abs(s))
plot(abs(n))


S = max(abs(s));
N = rms(n);

SNR = S/N
SNR_perRtScans = SNR/sqrt(2*ap.ns)

%% Plot T2D data

figure(1)
surf(echoVec/1000,1000*deltaVec',data)
shading flat
xlabel('T2 [ms]')
ylabel('delta [ms]')

figure(2)
surf(echoVec/1000,xD',dataY)
shading flat
xlabel('T2 [ms]')
set(gca,'defaulttextinterpreter','latex')
ylabel('$-\gamma^{2}G^{2}\delta^{2}(\Delta+\frac{\delta}{3})\times 10^{-9}$')
