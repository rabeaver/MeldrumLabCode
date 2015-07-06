clear
clc
close all

[spec,spec2,spec3] = readTecmag4d('C:\Users\tkmeldrum\Desktop\H2O_gradient_02.tnt');
nEchoes = 33;
nPts = 69;
nOmit = 5;
n2DPts = 11;

data = reshape(spec3,n2DPts,nEchoes,nPts);
data = data(:,:,1:end-nOmit);
dataRe = real(data);
dataIm = imag(data);
dataRe = reshape(dataRe,n2DPts,33*64);

point = 178;

dataDecay = dataRe(:,point);
dataNorm = dataDecay./dataDecay(1);
dataLog = log(dataNorm);
plot(dataLog)
% use points 1-6


gamma = 267.513e6; %rad s-1 T-1
D = 2.154e-9; %m2 s-1
delta = linspace(45e-6,375e-6,n2DPts);
DELTA = 2500e-6; %s

linPoints = -(gamma.*delta).^2.*D.*(DELTA + (2/3)*delta);

cftool(linPoints(1:6),dataLog(1:6));
