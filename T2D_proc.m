clear
clc
close all

%%

datadir = 'C:\CommonData\Membranes\PureWater\DELTAseries_Overnight_14July2016\';
datafile = 'SSET2Trad_membrane_PureWate__DELTA20000_14July2016_Overnight';


nPts = 56;                          % # of acqu points
omitPts = 0;                        % the number of points that are zeros from the spectrometer
nEchoes = 512;                      % Echoes
omitEchoes = 0;                     % numner of echoes to remove from data
tD = 2e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 200;                           % us
deltaMin = 20e-6;                  % s
deltaMax = 400e-6;                 % s
DELTA = 20e-3;                      % s
noisePoints = 5;                   % number of points for measuring noise
noiseNumber = 1;                    % scan number to use for determining SNR
G = 6.59;                           % T m-1, B0 field gradient


%%

gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1



[ap , spec, spec2,spec3,spec4] = readTecmag4d(strcat(datadir,datafile,'.tnt'));
T2Ddat = reshape(spec2, ap.td(2), nPts, nEchoes);
T2Ddat = T2Ddat(:,1:nPts-omitPts,omitEchoes+1:end);

deltaVec = linspace(deltaMin,deltaMax,ap.td(2));
BigDELTA = DELTA + deltaVec;
qIndex = 2*pi*gamma*1e6*G*deltaVec;
vIndex = (qIndex.^2.*(BigDELTA-deltaVec./3).*1e-9)';
% xD = -gammaRad^2*G^2.*deltaVec.^2.*(DELTA+2*deltaVec/3)*1e-9;

echoVec = tE*(omitEchoes+1):tE:(nEchoes*tE);

data = sum(real(T2Ddat),2);
% data = max(real(T2Ddat),[],2);
data = reshape(data,ap.td(2),(nEchoes-omitEchoes));
data = abs(data);

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
SNR_perRtScans = SNR/sqrt(ap.td(2)*ap.ns)

%% Plot T2D data

figure(1)
surf(echoVec/1000,vIndex',data)
shading flat
xlabel('T2 [ms]')
ylabel('delta [ms]')

t2axis = echoVec'*1e-6; %s

save(strcat(datadir,datafile, '.dat'), 'data', '-ascii')
save(strcat(datadir,datafile, '_T2axis.dat'), 't2axis', '-ascii')
save(strcat(datadir,datafile, '_vaxis.dat'), 'vIndex', '-ascii')

%UF Points [Min, Max; min(echoVec), max(echoVec), delta(eff)(min) [us], delta(eff)(max) [us], #echoes, #D points]
% sprintf('%f; %d %d %d; %.0f %.0f %.0f %.0f; %d %d',SNR, minind, maxind, firstinvertedind,  min(echoVec), max(echoVec), 1e6*min(t1), 1e6*max(t1), size(T1T2data,2), size(T1T2data,1))

% %%
% 
% fileID = fopen(strcat(datadir,'DataNotesAuto.txt'),'a');
% fprintf(fileID,'%s: %f; %.0f %.0f %.2f %.2f; %d %d\n',datafile, SNR, min(echoVec), max(echoVec), 1e6*min(deltaVec), 1e6*max(deltaVec), size(data,2), size(data,1));
% fclose(fileID);