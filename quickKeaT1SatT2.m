clear
clc
close all

% USER-DEFINED PARAMETERS
filename = 'data2D.2d';
T1axisFile = 'T1Axis.dat';
filedir = 'Z:\Data\TKM\WMO_2017\T1T2\P250.2015f\1\';

omitEchoes = 0; %front-end echoes to omit
% END USER-DEFINED PARAMETERS


fileloc = strcat(filedir,filename);
parloc  = strcat(filedir,'acqu.par');
T1axis = load(strcat(filedir,T1axisFile));

[ap,spec] = readKea4d(fileloc);
tE = readpar_Kea(parloc,'echoTime')*1e-6;
nrEchoes = ap.xDim;
nr2DPts = ap.yDim;

echoVec = (omitEchoes+1)*tE:tE:nrEchoes*tE;
spec2d = reshape(spec,nrEchoes,nr2DPts);
spec2d = spec2d(:,omitEchoes+1:end);
T1T2data = real(spec2d);

figure
surf(T1axis,echoVec(:,1:end)'*1000,T1T2data(:,1:end)); 
shading flat;
colormap('jet');
% shading interp;
colorbar 
xlabel('{\it t}_1 (ms)'); 
ylabel('{\it t}_2 (ms)');
title('T1-T2 data')


%% Save data, display ILT Data params
close all

T1T2data = T1T2data';

save(strcat(filedir,'T1T2data.out'), 'T1T2data', '-ascii');
save(strcat(datadir,datafile,'_T2axis.dat'),'T2axis','-ascii')
save(strcat(datadir,datafile,'_vaxis.dat'),'T1axis','-ascii')
