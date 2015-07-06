clear
close all
clc

%%

dir = 'C:\Users\tkmeldrum\Dropbox\Data\GdWater_T2D\';
cd(dir)

[aq,spec,spec2] = readTecmag4d('GdWater_2Oct2014.tnt');

spec3 = reshape(spec2',69,256,25);

sumData = sum(real(spec3),1);
sumData = reshape(sumData,256,25);
sumData = sumData';

tE = 150e-6;
nE = 256;
eV = (tE:tE:tE*nE)';

tau = linspace(0.1,0.036*25+0.1,25)';

save('data.dat','sumData','-ascii');
save('eV.dat','eV','-ascii');
save('tau.dat','tau','-ascii');