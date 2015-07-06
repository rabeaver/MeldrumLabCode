function [yphased,phaseangle] = autophase(ydata,phaseResolution)

phaseList = (-180:phaseResolution:180-phaseResolution)';
phaseListRad = phaseList*(2*pi/360);

yphased = zeros(size(ydata,1),size(ydata,2));
phaseangle = zeros(size(ydata,2),1);

for p = 1:size(ydata,2)
    cplxData = ydata(:,p);
    realData = real(ydata(:,p));
 
    cplxAngle = angle(cplxData);

    phasedDataList = zeros(size(cplxData,1),size(phaseList,1));
    newCplxAngleList = zeros(size(cplxData,1),1);


    for i = 1:1:size(phaseList,1)
        newCplxAngleList(:) = cplxAngle + phaseListRad(i);
        phasedDataList(:,i) = abs(cplxData).*(exp(1i*newCplxAngleList(:)));
    end

    phasedDataRealSum = sum(real(phasedDataList),1);
    [~,I] = max(phasedDataRealSum);
    
    phaseAngleChange = phaseListRad(I);
    
    yphased(:,p) = abs(cplxData).*exp(1i*(cplxAngle+phaseAngleChange));
    
    phaseangle(p) = phaseList(I);
end





