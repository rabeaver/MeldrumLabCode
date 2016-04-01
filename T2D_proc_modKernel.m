clear
clc
close all

%%
datadir = '/Users/tyler/Google Drive/Data2016/Tecmag/Acetone/';
datafile = 'AcetoneLarge_STE_31Mar2016_1';


nPts = 54;                          % # of acqu points
omitPts = 4;                        % the number of points that are zeros from the spectrometer
nEchoes = 64;                      % Echoes
omitEchoes = 2;                     % numner of echoes to remove from data
tD = 6e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 400;                           % us
deltaMin = 20e-6;                  % s
deltaMax = 600e-6;                 % s
DELTA = 0.5e-3;                      % s
noisePoints = 4;                   % number of points for measuring noise
noiseNumber = 1;                    % scan number to use for determining SNR
G = 6.59;                           % T m-1, B0 field gradient

%new code to account for fixed Td (total diffusion time)
Td = 1700e-6; %s

gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1


[ap , spec, spec2,spec3,spec4] = readTecmag4d(strcat(datadir,datafile,'.tnt'));
T2Ddat = reshape(spec2, ap.td(2), nPts, nEchoes);
T2Ddat = T2Ddat(:,1:nPts-omitPts,omitEchoes+1:end);

deltaVec = linspace(deltaMin,deltaMax,ap.td(2));
%new code to account for fixed Td (total diffusion time)
DELTAVec = linspace(Td-2*deltaMin,Td-2*deltaMax,ap.td(2));
DELTA = DELTAVec + deltaVec;

xD = gammaRad^2*G^2.*deltaVec.^2.*(DELTA-deltaVec/3);
% xD = gammaRad^2*G^2.*deltaVec.^2.*(DELTA+2*deltaVec/3);

echoVec = (tE*(omitEchoes+1):tE:(nEchoes*tE))*1e-6;

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
surf(echoVec,xD',data)
shading flat
xlabel('T2 [s]')
ylabel('delta [s]')

save(strcat(datadir,datafile, '.dat'), 'data', '-ascii')
% csvwrite(strcat(datadir,datafile, '.csv'), data)
save(strcat(datadir,datafile, '_T2dim.dat'), 'echoVec', '-ascii')
save(strcat(datadir,datafile, '_DDim.dat'), 'xD', '-ascii')

%UF Points [Min, Max; min(echoVec), max(echoVec), delta(eff)(min) [us], delta(eff)(max) [us], #echoes, #D points]
% sprintf('%f; %d %d %d; %.0f %.0f %.0f %.0f; %d %d',SNR, minind, maxind, firstinvertedind,  min(echoVec), max(echoVec), 1e6*min(t1), 1e6*max(t1), size(T1T2data,2), size(T1T2data,1))

fileID = fopen(strcat(datadir,'DataNotesAuto.txt'),'a');
fprintf(fileID,'%s: %f; %.0f %.0f %.2f %.2f; %d %d\n',datafile, SNR, min(echoVec), max(echoVec), 1e6*min(deltaVec), 1e6*max(deltaVec), size(data,2), size(data,1));
fclose(fileID);


