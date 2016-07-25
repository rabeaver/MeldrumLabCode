clear
clc
close all

%% Import Kea Data

filename = 'data2D_Cp';
fileext = '.2d';
filedir = 'C:\Users\tkmeldrum\Desktop\P250_2014e_SSET2Trad_6July2016\1\';

fileloc = strcat(filedir,filename,fileext);           % String together file name
parloc  = strcat(filedir,'acqu.par');         % String together Acquisition param location

[ap,spec] = readKea4d(fileloc);               % import Data and ap
tE = readpar_Kea(parloc,'echoTime')*1e-6;     % Set Echotime (s)
DELTA = readpar_Kea(parloc,'DELTA')*1e-3;     % Big DELTA (s)
deltamin = readpar_Kea(parloc,'tauMin')*1e-3; % (s)
deltamax = readpar_Kea(parloc,'tauMax')*1e-3; % (s)
nEchoes = ap.xDim;                            % Number of Echoes
nPtsD = ap.yDim;                              % Number of points
deltavec = linspace(deltamin,deltamax,nPtsD)';
gamma = 42.576;                               % MHz T-1
G = 23.87;                                    % T/m
tD = readpar_Kea(parloc,'dwellTime')*1e-6;    % dwell time (Tecmag shows correct dwell time for a complex point, no need to multiply by 2)

echoVec = tE:tE:nEchoes*tE;                   % Make Echovector
echoVec = echoVec';
spec2d = reshape(spec,nEchoes,nPtsD);         % Reshape Data
fitdata = sum(real(spec2d),1);                % Sum the data (Probably won't need for ILT)





%% Calculate other stuff

% Need to only use points where abs(fIndex)<=BWchirp/2
BigIndex = 1:nPtsD;

FOV = 1/(gamma*1e6*G*tD);
m_per_pt = FOV/nPtsD;
BigDELTA = DELTA + deltamax;

ptIndex = (1:nPtsD);
zBigIndex = FOV/2-BigIndex*m_per_pt;
fBigIndex = zBigIndex * gamma*1e6 * G;

qIndex = 2*pi*gamma*1e6*G*deltavec;
vIndex = qIndex.^2.*(BigDELTA-deltavec).*1e-9;

vIndex = rot90(vIndex,2);
vIndex = flipud(vIndex);
T2Ddat = (abs(spec2d))';



%% Save Data

save(strcat(filedir,filename,'_ILT','.dat'), 'T2Ddat', '-ascii')
save(strcat(filedir,filename, '_T2axis.dat'), 'echoVec', '-ascii')
save(strcat(filedir,filename, '_vaxis.dat'), 'vIndex', '-ascii')

