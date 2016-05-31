clear
clc
close all

filename = 'data.1d';
filedir = 'Z:\Data\TKM\GCI_Acrylate\17May2016\Rubber_CPMG\1\';
fileloc = strcat(filedir,filename);
parloc  = strcat(filedir,'acqu.par');

[ap,spec] = readKea4d(fileloc);


guess = [0.4 1 7]; %T2 in ms
[beta,R,J,CovB] = nlinfit(1e-3*spec(:,1),spec(:,2)./spec(1,2), @t2bifit_ampSumFixed, guess);
ci = nlparci(beta,R,'jacobian',J);

ypred = t2bifit_ampSumFixed(beta,1e-3*spec(:,1));

hh = figure(1);
hold on
scatter(1e-3*spec(:,1),spec(:,2)./spec(1,2));
plot(1e-3*spec(:,1),ypred,'-r');
xlabel('time/ms')
% pubgraph(hh, 16, 2, 'w')

sprintf('T2 = %f +/- %.1g ms.',beta(2), beta(2)-ci(2,1))