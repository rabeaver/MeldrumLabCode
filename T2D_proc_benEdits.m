clear
clc
close all

%%
datadir = 'C:\CommonData\EVOO\';
datafile = 'EVOOLarge_TraditSTE_100-2000de_15mDEL_16Mar2016';


nPts = 30;                          % # of acqu points
omitPts = 4;                        % the number of points that are zeros from the spectrometer
nEchoes = 512;                      % Echoes
omitEchoes = 2;                     % numner of echoes to remove from data
tD = 20e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 700;                           % us
deltaMin = 100e-6;                  % s
deltaMax = 2000e-6;                 % s
DELTA = 15e-3;                      % s
noisePoints = 1;                   % number of points for measuring noise
noiseNumber = 1;                    % scan number to use for determining SNR
G = 6.59;                           % T m-1, B0 field gradient
t2Sample = 85e-3;                   % s
t1Sample = 88e-3;                   % s



gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1


[ap , spec, spec2,spec3,spec4] = readTecmag4d(strcat(datadir,datafile,'.tnt'));
T2Ddat = reshape(spec2, ap.td(2), nPts, nEchoes);
T2Ddat = T2Ddat(:,1:nPts-omitPts,omitEchoes+1:end);


%%

deltaVec = linspace(deltaMin,deltaMax,ap.td(2));
xD = -gammaRad^2*G^2.*deltaVec.^2.*(DELTA+2*deltaVec/3)*1e-9;

echoVec = tE*(omitEchoes+1):tE:(nEchoes*tE);

data = sum(real(T2Ddat),2);
data = reshape(data,ap.td(2),(nEchoes-omitEchoes));

dataY = zeros(ap.td(2),nEchoes-omitEchoes);

for i = 1:ap.td(2)
    dataY(i,:) = log10(data(i,:)./data(i,1));
end

%% Other things
T2Ddata = sum(sum(T2Ddat,2),3); %crops data set according to above indices
T2Ddata = real(T2Ddata);

yD = log(T2Ddata./T2Ddata(1))';
xDif = -gammaRad^2*G^2.*deltaVec.^2.*(DELTA+2*deltaVec/3);
yCorr = yD +2*(deltaVec)./t2Sample +DELTA/t1Sample;

T2Dsize = size(T2Ddat,1); % cuts down delta points to math those selected for the indices,assuming that the 

% fit to STE diffusion equation
p = polyfit(xDif(1:T2Dsize),(yD(1,:)),1);
p2 = polyfit(xDif(1:T2Dsize),(yCorr(1,:)),1);

figure(9)
hold on
plot(xDif(1:T2Dsize),polyval(p,xDif(1:T2Dsize)),'k');
plot(xDif(1:T2Dsize),polyval(p2,xDif(1:T2Dsize)),'r');
scatter(xDif(1:T2Dsize),(yD(1,:)),'k') 
scatter(xDif(1:T2Dsize),(yCorr(1,:)),'r') 


legend('uncorrected','Corrected')
D = p(1)        % m^2 s^-1
Dcorr = p2(1)
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

save(strcat(datadir,datafile, '.dat'), 'data', '-ascii')

%UF Points [Min, Max; min(echoVec), max(echoVec), delta(eff)(min) [us], delta(eff)(max) [us], #echoes, #D points]
% sprintf('%f; %d %d %d; %.0f %.0f %.0f %.0f; %d %d',SNR, minind, maxind, firstinvertedind,  min(echoVec), max(echoVec), 1e6*min(t1), 1e6*max(t1), size(T1T2data,2), size(T1T2data,1))

fileID = fopen(strcat(datadir,'DataNotesAuto.txt'),'a');
fprintf(fileID,'%s: %f; %.0f %.0f %.2f %.2f; %d %d\n',datafile, SNR, min(echoVec), max(echoVec), 1e6*min(deltaVec), 1e6*max(deltaVec), size(data,2), size(data,1));
fclose(fileID);

