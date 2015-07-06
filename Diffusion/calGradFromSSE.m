close all
clc
clear

%%
% Enter the directory for the SSE data set on DI Water, as well as the
% temperature and the maximum data point to use.
dir = '/Users/tyler/Dropbox/Data/GradCal/1/';
T = 22.25; %degC
maxpt = 20;

%%
cd(dir)
data = load('data.dat');
data(:,1) = data(:,1); %time

gamma = 2.675222005e8; %s-1 T-1
gammaM = 42.57748; %MHz/T
D = selfDiffWater(T); %m^2/s
DEL = readpar_Kea('acqu.par','DELTA'); %read from acqu.par


tauAxis2 = -1e-9*(gamma^2*D).*(DEL.*data(:,1).^2+(2/3)*data(:,1).^3);

lnData = log(data(:,3)./data(1,3)); %ln signal intensity

scatter(tauAxis2(1:maxpt),lnData(1:maxpt))

% xaxis = -gamma^2*D*data(:,1).^2.*(DEL+(2/3)*data(:,1));

[slope,err] = polyfit(tauAxis2(1:maxpt),lnData(1:maxpt),1);
G = sqrt(slope(1))*gammaM;

res = lnData(1:maxpt)-slope(1)*tauAxis2(1:maxpt)-slope(2);
nSx2 = maxpt*sum(tauAxis2(1:maxpt).^2);
Sx2 = (sum(tauAxis2(1:maxpt)))^2;

S = sqrt(sum(res.^2)/(maxpt-2));
slopeErr = S*sqrt(maxpt/(nSx2-Sx2));

Gerr = 0.5*G^-0.5*slopeErr*gammaM;


% figure (1)
% hold on
% scatter(xaxis(1:maxpt),lnData(1:maxpt))

%%
GkHz = G/gammaM; %G in T/m

xaxisNew = -1e-9*gamma^2*G^2*data(:,1).^2.*(DEL+(2/3)*data(:,1));
figure(2)
plot(xaxisNew(1:maxpt),lnData(1:maxpt))

[slopeNew] = polyfit(xaxisNew(1:maxpt),lnData(1:maxpt),1);

D = 1e-9*slopeNew(1)
