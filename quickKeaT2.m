clear
clc
close all

filename = 'data.1d';
filedir = '/Volumes/ISC1026/Data/TKM/NGA_12May2016/Hamada/CDD_CPMG/1/';
fileloc = strcat(filedir,filename);
parloc  = strcat(filedir,'acqu.par');

[ap,spec] = readKea4d(fileloc);


guess = [0.4 1]; %T2 in ms
[beta,R,J,CovB] = nlinfit(1e-3*spec(:,1),spec(:,2), @t2monofit_simple, guess);
ci = nlparci(beta,R,'jacobian',J);

ypred = t2monofit_simple(beta,1e-3*spec(:,1));

hh = figure(1);
hold on
scatter(1e-3*spec(:,1),spec(:,2));
plot(1e-3*spec(:,1),ypred,'-r');
xlabel('time/ms')
% pubgraph(hh, 16, 2, 'w')

sprintf('T2 = %f +/- %.1g ms.',beta(2), beta(2)-ci(2,1))