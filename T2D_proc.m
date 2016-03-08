clear
clc
close all

%%
datadir = 'C:\CommonData\EVOO\';
datafile = 'EVOOLarge_TraditSTE_15mDEL_5Feb2016_result1';


nPts = 76;                          % # of acqu points
omitPts = 4;                        % the number of points that are zeros from the spectrometer
nEchoes = 256;                      % Echoes
omitEchoes = 4;                     % numner of echoes to remove from data
tD = 4e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 500;                           % us
deltaMin = 100e-6;                  % s
deltaMax = 1000e-6;                 % s
DELTA = 15e-3;                      % s
noisePoints = 4;                   % number of points for measuring noise
noiseNumber = 1;                    % scan number to use for determining SNR
G = 6.59;                           % T m-1, B0 field gradient



gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1


[ap , spec, spec2,spec3,spec4] = readTecmag4d(strcat(datadir,datafile,'.tnt'));
T2Ddat = reshape(spec2, ap.td(2), nPts, nEchoes);
T2Ddat = T2Ddat(:,1:nPts-omitPts,omitEchoes+1:end);

deltaVec = linspace(deltaMin,deltaMax,ap.td(2));
xD = -gammaRad^2*G^2.*deltaVec.^2.*(DELTA+deltaVec/3)*1e-9;

echoVec = tE*(omitEchoes+1):tE:(nEchoes*tE);

data = sum(real(T2Ddat),2);
data = reshape(data,ap.td(2),(nEchoes-omitEchoes));

dataY = zeros(ap.td(2),nEchoes-omitEchoes);

for i = 1:ap.td(2)
    dataY(i,:) = log10(data(i,:)./data(i,1));
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

save(strcat(datadir,datafile, '.dat'), 'data', '-ascii')

%UF Points [Min, Max; min(echoVec), max(echoVec), delta(eff)(min) [us], delta(eff)(max) [us], #echoes, #D points]
% sprintf('%f; %d %d %d; %.0f %.0f %.0f %.0f; %d %d',SNR, minind, maxind, firstinvertedind,  min(echoVec), max(echoVec), 1e6*min(t1), 1e6*max(t1), size(T1T2data,2), size(T1T2data,1))

fileID = fopen(strcat(datadir,'DataNotesAuto.txt'),'a');
fprintf(fileID,'%s: %f; %.0f %.0f %.2f %.2f; %d %d\n',datafile, SNR, min(echoVec), max(echoVec), 1e6*min(deltaVec), 1e6*max(deltaVec), size(data,2), size(data,1));
fclose(fileID);
