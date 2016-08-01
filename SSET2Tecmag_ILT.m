clear
clc
close all

%% User parameters
filename = '3per_CuH2O_SSET2Trad_19July2016_DELTAseries.tnt';
filedir = 'C:\CommonData\Membranes\CopperDELTAseriescheck\';
fileloc = strcat(filedir,filename);

[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);
tE = 200; %us
nEchoes = 512;
nPts = 56;
tD = 2e-6;
nPtsD = 21;
nPtsBlank = 0;
omitEchoes = 0;
DELTA = 3e-3; %s
deltamin = 20e-6; %s
deltamax = 400e-6; %s
deltavec = linspace(deltamin,deltamax,nPtsD)';
G = 6.59; %MHz m-1
gamma = 42.576;

data = reshape(spec2',nPts,nEchoes,nPtsD);
data = data(1:(nPts-nPtsBlank),omitEchoes+1:nEchoes,:);
dataInt = sum(data,1);
spec2d = reshape(dataInt,nEchoes,nPtsD);

echoVec = tE:tE:nEchoes*tE;                   % Make Echovector
echoVec = echoVec'./1e6;         % Reshape Data



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
T2Ddat = abs(spec2d)';



%% Save Data

save(strcat(filedir,filename,'_ILT','.dat'), 'T2Ddat', '-ascii')
save(strcat(filedir,filename, '_T2axis.dat'), 'echoVec', '-ascii')
save(strcat(filedir,filename, '_vaxis.dat'), 'vIndex', '-ascii')

