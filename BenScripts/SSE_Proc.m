clear
close all
clc
addpath(genpath('Z:\TKM\'));

%%Armory Processing
maindir = 'C:\Users\bmfortman\Documents\Data\ArmorySSE\';%remember to change both
cd(maindir);
adirnums = [2 3 5 6 7 8 9 10 11 12 14 15 16 17]; % for armory
atime = [1 10 12 14 30 33 35 37 39 42 44 46 49 51] % for armory change when adding more data

for j = 1:length(adirnums)
    dir = strcat(maindir,num2str(adirnums(j)),'\');
    cd(dir);
    
data1 = load('data.dat');
tau = ((data1((2 : 15),2))*10^9);
modsig = data1(:,3);
y1 = ((modsig(2 : 15))/(modsig(1)));
y = log(y1);
sig = data1(:,2);
figure(j)
hold on
plot(tau,y);

params(j).echoTime = readpar_Kea('acqu.par','echoTime');
params(j).pulseLength = readpar_Kea('acqu.par','pulseLength');
params(j).dwellTime = readpar_Kea('acqu.par','dwellTime');
params(j).nrEchoes = readpar_Kea('acqu.par','nrEchoes');
params(j).nrPnts = readpar_Kea('acqu.par','nrPnts');
params(j).nrScans = readpar_Kea('acqu.par','nrScans');
params(j).b1Freq = readpar_Kea('acqu.par','b1Freq');

%s0 = sig(1);
%y = log(sig(2:end)/s0);

%b = modsig(1);
%a1 = (modsig/(tau-b));
%a = sum(a1(:,15));
%plot(tau,a+b);

[p,S] = polyfit(tau,y,1);
fit= p(1)*tau+p(2); 
plot(tau,fit,'-r');
[fitci(j).ypred,fitci(j).delta] = polyval(p,tau,S) %this needs to be adjusted so that the % error in the slopes (as compared to real values can be determined)

sx = sum(tau);%tyler's error code
sy = sum(y);
sxx = sum(tau.^2);
sxy = sum(y.*tau);
syy = sum(y.^2);
N = length(tau);
D = N*sxx-sx^2;

a0(j) = (sxx*sy-sx*sxy)/D;
a1(j) = (N*sxy-sx*sy)/D;

stdErra(j).a0 = sqrt((sxx*(syy-a1(j)*sxy-a0(j)*sy)/((N-2)*D)));
stdErra(j).a1 = sqrt(N/sxx)*stdErra(j).a0;%Tyler's error code

xlabel('tau')
ylabel('modsig')
legend('data','line of best fit')
aslopes(j) = (p(1));
aslopesb(j) = (p(2));
clear('data1','tau','modsig','sig','p','y1','y')
end

%%Jamestown Processing
maindir = 'C:\Users\bmfortman\Documents\Data\JamestownSSE\';
cd(maindir);

jdirnums = [1 : 13]; %for jamestown
jtime = [1 10 12 14 16 30 35 39 42 44 46 49 51] %for jamestown
for j = 1:length(jdirnums)
    dir = strcat(maindir,num2str(jdirnums(j)),'\');
    cd(dir);
    
data1 = load('data.dat');
tau = ((data1((2 : 15),2))*10^9);
modsig = data1(:,3);
y1 = ((modsig(2 : 15))/(modsig(1)));
y = log(y1);
sig = data1(:,2);
figure(j)
hold on
plot(tau,y);

params(j).echoTime = readpar_Kea('acqu.par','echoTime');
params(j).pulseLength = readpar_Kea('acqu.par','pulseLength');
params(j).dwellTime = readpar_Kea('acqu.par','dwellTime');
params(j).nrEchoes = readpar_Kea('acqu.par','nrEchoes');
params(j).nrPnts = readpar_Kea('acqu.par','nrPnts');
params(j).nrScans = readpar_Kea('acqu.par','nrScans');
params(j).b1Freq = readpar_Kea('acqu.par','b1Freq');

[p,S] = polyfit(tau,y,1);
fit= p(1)*tau+p(2);
plot(tau,fit,'-r');
[fitci(j).ypred,fitci(j).delta] = polyval(p,tau,S)

sx = sum(tau);%tyler's error code
sy = sum(y);
sxx = sum(tau.^2);
sxy = sum(y.*tau);
syy = sum(y.^2);
N = length(tau);
D = N*sxx-sx^2;

a0(j) = (sxx*sy-sx*sxy)/D;
a1(j) = (N*sxy-sx*sy)/D;

stdErrj(j).a0 = sqrt((sxx*(syy-a1(j)*sxy-a0(j)*sy)/((N-2)*D)));
stdErrj(j).a1 = sqrt(N/sxx)*stdErrj(j).a0;%Tyler's error code

xlabel('tau')
ylabel('modsig')
legend('data','line of best fit')
jslopes(j) = (p(1));
jslopesb(j) = (p(2));
clear('data1','tau','modsig','sig','p','y1','y')
end


%%final slope comparisons
close all
figure (14)
hold on
atime = atime';
jtime = jtime';
for j = 1:length(jdirnums)
    errSlpej(j) = stdErrj(j).a1;
end
for j = 1:length(adirnums)
    errSlpea(j) = stdErra(j).a1;
end
errorbar(jtime,jslopes,errSlpej,'-k')
errorbar(atime,aslopes,errSlpea,'-r')
legend('jamestown','armory')
ylabel('Diffusion of water')
xlabel('Days')
