clc
close all
clear

%%
% Info for Jamestown Mortar sample measured on 17 July 2014 by TKM
filename = '/Users/tyler/Dropbox/Data/MortarDiffusion/JamestownMortar_17July2014.tnt';

[ap,spec1,spec2] = readTecmag(filename);
nPts = 69;
oPts = 5;
td = ap.dw(1); %s
nScans = 256;
n2D = ap.td(2);
tE = 200e-6; %s
nEchoes = 256;
Td = 1500e-6; %s
delta = 1e-6*(87.5:12.5:700)';
DELTA = 1e-6*(1325:-25:100)';
gamma = 2.675222e8; % s-1 T-1
g = 280; %kHz mm-1
G = g/42.57748; %T m-1

spec3 = reshape(spec2',nPts,nEchoes,n2D);
spec3 = spec3(1:nPts-oPts,:,:);

dataInt = sum(spec3,1);
dataInt = reshape(dataInt,nEchoes,n2D)';
echoVec = (tE:tE:tE*nEchoes)';

q2 = (delta*gamma*G).^2;

save('echoVec.txt','echoVec','-ascii');
save('delta.txt','delta','-ascii');
save('q2.txt','q2','-ascii');
save('dataInt.txt','dataInt','-ascii');


