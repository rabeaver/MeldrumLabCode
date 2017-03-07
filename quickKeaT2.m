clear
clc
close all

% USER-DEFINED PARAMETERS
filename = 'data.2d';
filedir = 'Z:\Data\TKM\WMO_2017\CPMG\P250.2015f\5\';

omitEchoes = 0; %front-end echoes to omit
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


guess = [0 1 0.0002]; %T2 in s
[beta,R,J,CovB] = nlinfit(echoVec,fitdata./fitdata(1), @t2monofit, guess);
ci = nlparci(beta,R,'jacobian',J);

ypred = t2monofit(beta,echoVec);

hh = figure(1);
hold on
scatter(echoVec,fitdata./fitdata(1));
plot(echoVec,ypred,'-r');
xlabel('time/ms')
% pubgraph(hh, 16, 2, 'w')

sprintf('T2 = %f +/- %.1g ms.',1000*beta(3), 1000*beta(3)-1000*ci(3,1))