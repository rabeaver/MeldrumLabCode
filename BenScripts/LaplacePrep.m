%% starting the inverse laplace stuff
% close all
clear all

cd('C:\Users\bmfortman\Documents\Data\MortarDiffusion\')
[aq,spec,spec2]=readTecmag('Sample8_2D16_n64_16June2014.tnt'); % reads in data

spec3 = reshape(spec2',69,32,16); % reshapes to params, points, number echoes, tau points
spec3 = spec3(1:64,:,:); % cuts out blank points
intdata = sum(spec3,1); % integrates echoes
intdata = abs(intdata); % takes abs value of those integrals
intdata = reshape(intdata,32,16); % reshapes to a matrix; nrEchos, tau points
tEcho = (1:32)*150e-6; % creates time axis
tEcho = tEcho';
save('intdata.txt','intdata','-ascii'); % saves in a text file for use with 2dLaplace
save('tEcho.txt','tEcho','-ascii');% saves in a text file for use with 2dLaplace

% these tD files are taken from the excel sheet

tD1 = [9.000
37.000
65.000
93.000
121.000
149.000
177.000
205.000
233.000
261.000
289.000
317.000
345.000
373.000
401.000
429.000];

tD2 = [964.000
908.000
852.000
796.000
740.000
684.000
628.000
572.000
516.000
460.000
404.000
348.000
292.000
236.000
180.000
124.000];

tD3 = [87.000
115.000
143.000
171.000
199.000
227.000
255.000
283.000
311.000
339.000
367.000
395.000
423.000
451.000
479.000
507.000];

tD1 = tD1/1000;
tD2 = tD2/1000;
tD3 = tD3/1000;

save('tD1.txt','tD1','-ascii');
save('tD2.txt','tD2','-ascii');
save('tD3.txt','tD3','-ascii');

TwoDLaplaceInverse