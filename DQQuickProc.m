clear
clc
close all

%%
info(1).dir = '/Users/tyler/Dropbox/Data/DQ/DQ_RubberBand/22_expTime/';
info(2).dir = '/Users/tyler/Dropbox/Data/DQ/DQ_RubberBand/23_expTimerec0/';
paramsFile  = 'acqu.par';
reFile      = 'dataRe.dat';
imFile      = 'dataIm.dat';
dataFile    = 'data.dat';
% autophase   = 1;
nrEndPnts   =9
%% DQ
for k=1:2;
cd(info(k).dir)

dataRe = load(reFile);
dataIm = load(imFile);
datadat= load(dataFile);
tauVector = datadat(:,1);
dataCp = complex(dataRe,dataIm);
% 
% if autophase == 1
%     dataCp = autophase(dataCp',1);
% else
% dataCp = dataCp';
% end

params.acqTime  = readpar_Kea(paramsFile,'acqTime');
params.bandwidth  = readpar_Kea(paramsFile,'bandwidth');
params.nrScans  = readpar_Kea(paramsFile,'nrScans');
params.rxPhase  = readpar_Kea(paramsFile,'rxPhase');
params.rxGain  = readpar_Kea(paramsFile,'rxGain');
params.nrPts  = readpar_Kea(paramsFile,'nrPnts');
params.repTime  = readpar_Kea(paramsFile,'repTime');
params.b1Freq  = readpar_Kea(paramsFile,'b1Freq');
params.nrEchoes  = readpar_Kea(paramsFile,'nrEchoes');
params.echoTime  = readpar_Kea(paramsFile,'echoTime');
params.nrPntsTau  = readpar_Kea(paramsFile,'nrPntsTau');
params.nrPnts  = readpar_Kea(paramsFile,'nrPnts');
params.dwellTime  = readpar_Kea(paramsFile,'dwellTime');

echoVector = linspace(params.dwellTime,params.dwellTime*params.nrPnts,params.nrPnts);

for i = 1:params.nrPntsTau
    stack = reshape(dataCp(i,:),size(dataCp,2)/params.nrEchoes,params.nrEchoes);
    stacksum(:,i) = sum(stack,2);
end

output(k).dataRe = real(stacksum);
output(k).dataIm = imag(stacksum);

output(k).dataInt = sum(output(k).dataRe,1);

end


normFact = max(output(2).dataInt);
output(1).dataInt = output(1).dataInt/normFact;
output(2).dataInt = output(2).dataInt/normFact;
sumdata = output(1).dataInt + output(2).dataInt;

figure(1)
hold on
plot(tauVector,output(1).dataInt,'-k')
plot(tauVector,output(2).dataInt,'-b')
plot(tauVector,sumdata,'-r')


normdata = output(1).dataInt./sumdata;

A = ones(nrEndPnts,2);
A(:,1) = tauVector(1:nrEndPnts);
b = (output(2).dataInt(1:nrEndPnts))';
solution = A\b;

% output(2).subdata = output(2).dataInt - (solution(2)+solution(1)*tauVector');


% 
% output(2).subdata = output(2).subdata/normFact;

figure(4)
hold on
plot(tauVector,output(1).dataInt,'-ok')
plot(tauVector,output(2).dataInt,'-or')
% plot(tauVector,output(2).subdata,'-ob')

figure(5)
hold on
scatter(tauVector,output(2).dataInt)
plot(tauVector,(solution(2)+solution(1)*tauVector))
% plot(tauVector,(output(2).subdata),'or')

figure(6)
plot(tauVector,normdata)
% dataAb = abs(stacksum);
 %%
% figure(1)
% surf(tauVector,echoVector,dataRe)
% shading interp
% 
% figure(2)
% surf(tauVector,echoVector,dataIm)
% shading interp
% % 
% figure(3)
% surf(tauVector,echoVector,dataAb)
% shading flat
%%



%%

builduplim = 64;
figure(5)
plot(tauVector(1:builduplim),dataInt(1:builduplim))
%%
% FWHM = 7; 
% weight = exp((-(echoVector-params.acqTime*500).^2)/(2*FWHM/(2*sqrt(2*log(2)))));
% 
% 
% for i = 1:size(dataRe,2)
%     dataRe_w(:,i) = dataRe(:,i).*weight';
% end

%%
% figure(5)
% surf(tauVector,echoVector,dataRe_w)
% shading interp
% % zlim([-0.5 0.5])
% xlabel('echo Time')
% ylabel('tau time')
% view(90,90)
% 
% dataInt2 = sum(dataRe_w,1);
% figure(4)
% hold on
% plot(tauVector,dataInt2./max(dataInt2),'-r')
%%

figure(3)
surf(tauVector,echoVector,dataIm)
shading interp
zlim([-0.5 0.5])
% view(-90,0)

figure(4)
surf(tauVector,echoVector,dataAb)
shading interp
zlim([-0.5 0.5])
view(90,90)


%%
figure(5)
sumdata = sum(dataRe,2);
plot(tauVector,sumdata);
% % 
figure(5)
plot(tauVector,p)


