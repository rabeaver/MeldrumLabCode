clear
clc
close all

% USER-DEFINED PARAMETERS
filename = 'data.2d';
filedir = 'Z:\JNK\PM5\UFT2Ddata\Paint_Nick\LinseedOil_T2_28June2016\1\';

omitEchoes = 2; %front-end echoes to omit
% END USER-DEFINED PARAMETERS


fileloc = strcat(filedir,filename);
parloc  = strcat(filedir,'acqu.par');

[ap,spec] = readKea4d(fileloc);
tE = readpar_Kea(parloc,'echoTime')*1e-6;
nrEchoes = ap.yDim;
nrPts = ap.xDim;

echoVec = (omitEchoes+1)*tE:tE:nrEchoes*tE;
spec2d = reshape(spec,nrPts,nrEchoes);
spec2d = spec2d(:,omitEchoes+1:end);
fitdata = sum(real(spec2d),1);


guess = [0.4 1 7]; %T2 in ms
[beta,R,J,CovB] = nlinfit(echoVec,fitdata./fitdata(1), @t2bifit_ampSumFixed, guess);
ci = nlparci(beta,R,'jacobian',J);

ypred = t2bifit_ampSumFixed(beta,echoVec);

hh = figure(1);
hold on
scatter(echoVec,fitdata./fitdata(1));
plot(echoVec,ypred,'-r');
xlabel('time/ms')
% pubgraph(hh, 16, 2, 'w')

sprintf('T2 = %f +/- %.1g ms.',beta(2), beta(2)-ci(2,1))