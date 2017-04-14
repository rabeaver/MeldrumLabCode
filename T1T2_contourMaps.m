clear
clc
close all

%give file dir, file name (the *.out file from Prospa export2d), T1 and T2 limits (should be the same)
datadirM = '/Users/tyler/Dropbox/Data/WMO_2017/T1T2/M205.2015e/1/';
datafileM = '2DILT.out';
datadirP = '/Users/tyler/Dropbox/Data/WMO_2017/T1T2/P205.2015e/1/';
datafileP = '2DILT.out';
outdir = '/Users/tyler/Dropbox/Data/WMO_2017/T1T2/ComparisonFigs/';
outfileprefix = '205comp';
graphtitle = '205';
T2lims = [1e-4 1e-2];
T1lims = [1e-3 1e0];

%load the data and remove 0 values (replace with NaN)
dataM = load(strcat(datadirM,datafileM));
dataM = interp2(dataM,2);
dataM(dataM == 0) = NaN; 
nPts = size(dataM,1);

%load the data and remove 0 values (replace with NaN)
dataP = load(strcat(datadirP,datafileP));
dataP = interp2(dataP,2);
dataP(dataP == 0) = NaN; 

%this calculates the axes based on the limits above
T1axis = logspace(log10(T1lims(1)),log10(T1lims(2)),nPts);
T2axis = logspace(log10(T2lims(1)),log10(T2lims(2)),nPts);

save(strcat(outdir,outfileprefix,'.mat'))

%%plot the T1T2 data
hh = figure(1);
hold on
[CM,hm] = contour(T2axis,T1axis,dataM,7,'Color','k');
[CP,hp] = contour(T2axis,T1axis,dataP,7,'Color','r');
legend('M205','P205')
title(graphtitle);
set(gca,'XScale','log','YScale','log','XTick', [1e-4; 1e-3; 1e-2; 1e-1; 1e0; 1e1])
xlim(T2lims)
ylim(T1lims)
ylabel('\itT\rm_1 [s]')
xlabel('\itT\rm_2 [s]')
view([0,90])
pubgraph(hh,16,2,'w','Arial');
print(strcat(outdir,outfileprefix,'.jpg'),'-djpeg')
savefig(strcat(outdir,outfileprefix,'.fig'))



