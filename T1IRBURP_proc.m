clear
clc
close all  

%%

% Input filename, - .tnt
filename = 'T1IRBURP_degassedEthyleneGlycolMolecularSieves_19March2016_1024scans';
filedir = 'C:\CommonData\EthyleneGlycolMolecularSieves\';

fileloc = strcat(filedir,filename,'.tnt');

% Read file
[ap,spec,spec2,spec3,spec4] = readTecmag4d(fileloc);

% Input experiment parameters

tEcho = 700; %us
nEchoes = 512;
nPts = 76;
nPtsBlank = 1;
omitEchoes = 0; 
nT1Pts = 11;
T1min = .100; %ms
T1max = 2000; %ms
noisePoints = 1; %number of points to use for noise at beginning and end of each acqu period
noiseNumber = nT1Pts; %T1 point to use for SNR calc

echoVector = ((1+omitEchoes)*tEcho:tEcho:nEchoes*tEcho); % T2 vector


% Specify lin or log spaced points
linORlog = 0; % 0 for linearly space and 1 for log spaced

% Make T1vector
if linORlog == 0
    T1vector = linspace((T1min),(T1max),nT1Pts); % Linspace T1 points
else
    T1vector = logspace(log10(T1min),log10(T1max),nT1Pts); % Logspace T1sat
end
% %% SNR calc
% 
% fileloc = strcat(filedir,filename,'.tnt');
% 
% % Read file
% [ap,specN,spec,spec3,spec4] = readTecmag4d(fileloc);
% 
% 
% [~,Spoint] = max(abs(real(spec2(21,:))));
% Spoint = Spoint + 16*nPts;
% S = (real(spec2(nT1Pts,Spoint-nPts/2:Spoint+nPts/2)));
% N = (imag(spec2(nT1Pts,Spoint-nPts/2:Spoint+nPts/2)));
% % N = (real(specN(Spoint-nPts/2:Spoint+nPts/2)))';
% 
% SNR = snr(S,N)

%% SNR calc
data = reshape(spec2,nT1Pts,nPts,nEchoes);
data = data(:,1:nPts-nPtsBlank,omitEchoes+1:end);
n1 = data(noiseNumber,1:noisePoints,:);
n2 = data(noiseNumber,nPts-nPtsBlank-noisePoints:end,:);
ndata = cat(2,n1,n2);
ndata = reshape(ndata,1,(2*noisePoints+1)*(nEchoes-omitEchoes));
sdata = data(noiseNumber,:,:);
sdata = reshape(sdata,1,(nPts-nPtsBlank)*(nEchoes-omitEchoes));

% data = reshape(data,nT1Pts,(nPts-nPtsBlank)*nEchoes);

%
S = max(abs(sdata(end,:)));
N = rms(ndata(end,:));

% SNR = snr(S,N)
SNR = S/N
SNR_perRtScan = SNR/sqrt(nT1Pts*ap.ns)
% 
figure(1)
hold on
plot(abs(sdata))
plot(abs(ndata))
hold off

%% Make 2D data set for T1IRT2 ILT

data = reshape(spec2',nPts,nEchoes,nT1Pts);
data = data(1:(nPts-nPtsBlank),(1+omitEchoes):end,:);
data2d = sum(real(data),1);
data2d = reshape(data2d,nEchoes-omitEchoes,nT1Pts);
data2d = data2d';

% Plot of data
figure(2)
surf(echoVector,T1vector,data2d); shading flat

% Save data in specified directory with the same filename and ".dat"
% extension
save(strcat(filedir,filename,'.dat'), 'data2d', '-ascii')

 %% 1D Fits

%T1 (A*(1-2*exp(-x/T1))
cftool(T1vector, data2d(:,1)./abs(max(data2d(:,1))))


%%
%T2 (A*exp(-x/T2))
cftool(echoVector/10^3, data2d(end,:)'./max(data2d(end,:)))
