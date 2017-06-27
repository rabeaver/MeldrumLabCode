clear
clc
close all

%%
datadir = 'C:\CommonDataTKM\pH2\MeOH\';
datafile = 'CPMG_1_16June2017';


nPts = 78;                          % # of acqu points
omitPts = 0;                        % the number of points that are zeros from the spectrometer
nEchoes = 512;                      % Echoes
omitEchoes = 0;                     % numner of echoes to remove from data
tD = 1e-6;                          % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)
tE = 250;                           % us
deltaMin = 40e-6;                  % s
deltaMax = 2659e-6;                 % s
lin = 1;                            % 1 if delta is linearly spaced, 0 if log spaced
delta = 4e-3;                     % s
DELTA = 25e-3;                      % s
noisePoints = 2;                   % number of points for measuring noise
noiseNumber = 1;                    % scan number to use for determining SNR
G = 6.59;                           % T m-1, B0 field gradient


%%

gamma = 42.576;                     % MHz T-1
gammaRad = gamma*2*pi*1e6;          % rad s-1 T-1



[ap , spec, spec2,spec3,spec4] = readTecmag4d(strcat(datadir,datafile,'.tnt'));
T2Ddat = reshape(spec2, ap.td(2), nPts, nEchoes);
T2Ddat = T2Ddat(:,1:nPts-omitPts,omitEchoes+1:end);


if lin==1;
    deltaVec = linspace(deltaMin,deltaMax,ap.td(2));
else
    deltaVec = logspace(log10(deltaMin),log10(deltaMax),ap.td(2));
end

BigDELTA = DELTA + deltaVec;
qIndex = 2*pi*gamma*1e6*G*deltaVec;
vIndex = (qIndex.^2.*(BigDELTA-deltaVec./3).*1e-9)';
vIndex_RC = (1/6)*(gammaRad*G)^2.*(delta^3 - 3*delta^2*deltaVec/2 + 3*(deltaVec/2).^2*delta + 6*(deltaVec/2).^2*DELTA + 3*(deltaVec/2).^3)*1e-9;
vIndex_RC = vIndex_RC';
% xD = -gammaRad^2*G^2.*deltaVec.^2.*(DELTA+2*deltaVec/3)*1e-9;

echoVec = tE*(omitEchoes+1):tE:(nEchoes*tE);

data = sum(real(T2Ddat),2);
% data = max(real(T2Ddat),[],2);
data = reshape(data,ap.td(2),(nEchoes-omitEchoes));
data = real(data);

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

figure(1)
hold on
plot(abs(s))
plot(abs(n))



S = max(abs(s));
N = rms(n);

SNR = S/N
SNR_perRtScans = SNR/sqrt(ap.td(2)*ap.ns)

%% Plot T2D data

figure(2)
hold on
surf(echoVec/1000,deltaVec*1000,data./max(max(data)))
shading flat
xlabel('T2 [ms]')
ylabel('delta [ms]')
if lin==0
    set(gca,'YScale','log')
end

figure(3)
hold on
surf(echoVec/1000,vIndex_RC,data./max(max(data)))
shading flat
xlabel('T2 [ms]')
ylabel('v [s m-2]')
if lin==0
    set(gca,'YScale','log')
end

t2axis = echoVec'*1e-6; %s

save(strcat(datadir,datafile, '.dat'), 'data', '-ascii')
save(strcat(datadir,datafile, '_T2axis.dat'), 't2axis', '-ascii')
save(strcat(datadir,datafile, '_vaxis.dat'), 'vIndex', '-ascii')
save(strcat(datadir,datafile, '_vaxisRC.dat'), 'vIndex_RC', '-ascii')

%UF Points [Min, Max; min(echoVec), max(echoVec), delta(eff)(min) [us], delta(eff)(max) [us], #echoes, #D points]
% sprintf('%f; %d %d %d; %.0f %.0f %.0f %.0f; %d %d',SNR, minind, maxind, firstinvertedind,  min(echoVec), max(echoVec), 1e6*min(t1), 1e6*max(t1), size(T1T2data,2), size(T1T2data,1))

% %%
% 
% fileID = fopen(strcat(datadir,'DataNotesAuto.txt'),'a');
% fprintf(fileID,'%s: %f; %.0f %.0f %.2f %.2f; %d %d\n',datafile, SNR, min(echoVec), max(echoVec), 1e6*min(deltaVec), 1e6*max(deltaVec), size(data,2), size(data,1));
% fclose(fileID);