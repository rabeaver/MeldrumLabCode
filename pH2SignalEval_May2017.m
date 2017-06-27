clear
clc
close all

%%
filedir = 'C:\CommonData\TKM\pH2\';
filenameHOT = 'HOT_pH2_SignalTest_10May2017.tnt';
filenameCOLD = 'COLD_pH2_SignalTest_10May2017.tnt';

[apHOT,specHOT] = readTecmag4d(strcat(filedir,filenameHOT));
tHOT = apHOT.dw(1)*(1:1:apHOT.ctd);

[apCOLD,specCOLD] = readTecmag4d(strcat(filedir,filenameCOLD));
tCOLD = apCOLD.dw(1)*(1:1:apCOLD.ctd);

hh = figure(1);
hold on;
plot(tHOT,real(specHOT));
plot(tHOT,imag(specHOT));
plot(tCOLD,real(specCOLD));
plot(tCOLD,imag(specCOLD));
xlabel('acq time [s]');
ylabel('S [arb]')
legend('HOT real','HOT imag','COLD real','COLD imag')
pubgraph(hh,14,0.5,'w','Arial')

ii = figure(2);
hold on;
plot(tHOT,abs(specHOT),'-r');
plot(tCOLD,abs(specCOLD),'-b');
xlabel('acq time [s]');
ylabel('abs S [arb]')
legend('HOT','COLD')
pubgraph(ii,14,0.5,'w','Arial')