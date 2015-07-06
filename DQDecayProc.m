close all
clc
clear

%%
info(1).dir = '/Users/tyler/Dropbox/Data/DQ/DQ_Buildup/20/';
info(2).dir = '/Users/tyler/Dropbox/Data/DQ/DQ_Decay/20/';
nrIntEchoes = 32; %number of echoes to integrate together, between 1-total number of echoes
gFilter = 0; %apply Gaussian filter to each echo 1=yes, 0=no.

for k = 1:size(info,2);
cd(info(k).dir)

data(k).dat = load('data.dat');
data(k).dataRe = load('dataRe.dat');
data(k).dataIm = load('dataIm.dat');
data(k).dataCp = complex(data(k).dataRe,data(k).dataIm);
% data(k).dataRe_p = real(autophase(data(k).dataCp,1));
data(k).dataRe_p = data(k).dataRe;


params(k).echoTime = readpar_Kea('acqu.par','echoTime');
params(k).acqTime = readpar_Kea('acqu.par','acqTime');
params(k).evoTime = readpar_Kea('acqu.par','evoTime');
params(k).nrEchoes = readpar_Kea('acqu.par','nrEchoes');
params(k).nrPnts = readpar_Kea('acqu.par','nrPnts');
params(k).dwellTime = readpar_Kea('acqu.par','dwellTime');
params(k).nrPntsTau = readpar_Kea('acqu.par','nrPntsTau');
params(k).nrScans = readpar_Kea('acqu.par','nrScans');
params(k).pulseLength = readpar_Kea('acqu.par','pulseLength');
params(k).zFilt = readpar_Kea('acqu.par','zFilt');

data(k).tauVector = data(k).dat(:,1);

filter(k).x = params(k).dwellTime:params(k).dwellTime:params(k).nrPnts*params(k).dwellTime;
filter(k).c = params(k).pulseLength/(2*sqrt(2*log(2)));
filter(k).G = exp(-((filter(k).x-(params(k).dwellTime*params(k).nrPnts/2)).^2)/(2*filter(k).c^2));
filter(k).Gvec = repmat(filter(k).G,params(k).nrEchoes);
filter(k).Gvec = filter(k).Gvec(1:params(k).nrPntsTau,:);



if gFilter == 1
    data(k).filtData = data(k).dataRe_p.*filter(k).Gvec;
    data(k).dataInt = sum(data(k).filtData(:,1:nrIntEchoes*params(k).nrPnts),2);
elseif gFilter == 0
    data(k).filtData = data(k).dataRe_p;
    data(k).dataInt = sum(data(k).filtData(:,1:nrIntEchoes*params(k).nrPnts),2);
end

data(k).logDataInt = log10(data(k).dataInt);


end




figure(1)
hold on
plot(data(1).tauVector,data(1).dataInt,'-r')
plot(data(2).tauVector,data(2).dataInt,'-b')
%%
j = 1; %index for buildup 0202 curve
k = 2; %index for decay 0000 curve
lastEarlyPt = 58;
firstEndPt = 58;

[beta.early, S.early] = polyfit(data(k).tauVector(1:lastEarlyPt),data(k).logDataInt(1:lastEarlyPt),1);
[beta.late, S.late] = polyfit(data(k).tauVector(firstEndPt:end),data(k).logDataInt(firstEndPt:end),1);
lineData.early = polyval(beta.early,data(k).tauVector);
lineData.late = polyval(beta.late,data(k).tauVector);
res.early = data(k).logDataInt - lineData.early;
res.late = data(k).logDataInt - lineData.late;
data(k).subIntData = 10.^(res.late);

figure(2)
subplot(2,1,1)
hold on
scatter(data(k).tauVector,data(k).logDataInt)
plot(data(k).tauVector,lineData.early,'-b')
plot(data(k).tauVector,lineData.late,'-r')
ylabel('log integrated echo data')
subplot(2,1,2)
hold on
plot(data(k).tauVector,res.early,'ob')
plot(data(k).tauVector,res.late,'or')
line([0 max(data(k).tauVector)] ,[0 0],'Color',[0 0 0]);
ylim([-0.1 0.1])
ylabel('residual log integrated echo data from linear fits')
xlabel('tau time (ms)')

figure(3)
hold on
plot(data(j).tauVector,data(j).dataInt,'-r')
plot(data(k).tauVector,data(k).subIntData,'-b')
ylabel('subtracted integrated echo data')
xlabel('tau time (ms)')

normData = data(j).dataInt./(data(j).dataInt + data(k).dataInt);
figure(4)
plot(data(j).tauVector,normData)
ylim([-0.05 0.75])

%%

n=[1 32]; %which echoes to show
p = 23; %which tau point
figure(7)
plot(data(1).filtData(p,params(1).nrPnts*(n(1)-1)+1:params(1).nrPnts*n(2)))

%%
figure(5)
suptitle('Decay 0000')
for i=1:params(k).nrPntsTau;
    subplot(ceil(sqrt(params(k).nrPntsTau)),ceil(sqrt(params(k).nrPntsTau)),i)
    hold on
%     plot(data(j).filtData(i,:),'-r')
    plot(data(k).filtData(i,:),'-b')
    ylim([-200 400])
    xlim([0 params(k).nrPnts*params(k).nrEchoes])
end


figure(6)
suptitle('Buildup 0202')
for i=1:params(k).nrPntsTau;
    subplot(ceil(sqrt(params(k).nrPntsTau)),ceil(sqrt(params(k).nrPntsTau)),i)
    hold on
    plot(data(j).filtData(i,:),'-r')
%     plot(data(k).filtData(i,:),'-b')
    ylim([-40 40])
    xlim([0 params(k).nrPnts*params(k).nrEchoes])
end


%%
for j = 1:params.nrPntsTau
figure(1)
    for i = 1:params.nrEchoes
        subplot(round(sqrt(params.nrEchoes)),round(sqrt(params.nrEchoes)),i)
        plot(filtData(j,32*i-31:32*i))
        xlim([0 params.nrPnts*params.dwellTime])
        ylim([-20 20])
    end
end

%%

dataInt = sum(filtData,2);
figure(2)
plot(tauVector,dataInt)
