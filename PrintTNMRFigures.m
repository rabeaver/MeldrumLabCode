clear
clc

filename = ('PM25_CPMG_Rubber_powerModulation_5June2014_4us.tnt');
info = 'PM254us';

[aq,spec,spec2] = readTecmag(filename);

h = figure(1);
pcolor(abs(spec2));
shading flat
colormap gray
view([0 90])
xlim([ 0 2048])
ylim([0 50])
% set(gca,'YDir','reverse')
title(info)
xlabel('time [us]')
ylabel('power')
caxis(1e3*[0.2 .7])

print(h, '-djpeg', '-r300', strcat(info,'.jpg'))