close all
clear
clc

alldata = dlmread('T2data-1850.dat');

omitFirst = 1;

realData = alldata(omitFirst + 1:end,2);
imagData = alldata(omitFirst + 1:end,3);
timeaxis= alldata(omitFirst + 1:end,1);

cplxData = complex(realData,imagData);
cplxAngle = angle(cplxData);


%%

phaseList = (-180:0.1:179)';
phaseListRad = phaseList*(2*pi/360);

phase = 0;
phaseRad = phase*(2*pi/360);

newCplxAngle = cplxAngle + phaseRad;
phasedData = realData.*(exp(1i.*newCplxAngle));

phasedDataList = zeros(size(cplxData,1),size(phaseList,1));
newCplxAngleList = zeros(size(cplxData,1),size(phaseList,1));


for i = 1:1:size(phaseList,1)
    newCplxAngleList(:,i) = cplxAngle + phaseListRad(i);
    phasedDataList(:,i) = realData.*(exp(1i.*newCplxAngleList(:,i)));
end

% phasedData = realData.*(exp(1i.*newCplxAngle));


realPhased = real(phasedData);
imagPhased = imag(phasedData);

phasedDataRealSum = sum(real(phasedDataList),1);
[C,I] = max(phasedDataRealSum);

realBestPhased = real(phasedDataList(:,I));
imagBestPhased = imag(phasedDataList(:,I));

bestAngle = phaseList(I)

figure(1)
clf
hold on
scatter(timeaxis,realBestPhased,'ob')
scatter(timeaxis,imagBestPhased,'or')

% figure
% scatter(timeaxis,cplxAngle)