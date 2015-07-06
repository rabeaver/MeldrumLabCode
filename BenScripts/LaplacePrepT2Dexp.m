%% starting the inverse laplace stuff
close all
clear all

cd('C:\Users\bmfortman\Documents\Data\MortarDrying\')
[aq,spec,spec2]=readTecmag4d('Jamestown_69pts_256ec_150tE_T2D_1000Td_9-27-14.tnt'); % reads in data

nrEchoes = 256;
%parameters for the gaussian filter
pulseLength = 6;
acqTime = 69; % acqPts * dwellTime
dwellTime = 1; % in us
nrPts = acqTime; 




%% Gaussian filter


% for Gaussian
gaussianLb = 1;
t = dwellTime:dwellTime:(dwellTime*length(spec2)/nrEchoes);
% Creates Gaussian Matrices??
% g_x = exp(-(t-(max(t)/2)).^2/(2*(1/lb)^2)); % Supresses noise
c = pulseLength/(2*sqrt(2*log(2))); %parameters and fit for gaussian filter
g_x = exp(-((t-acqTime/2).^2)/(2*c^2));

g_x = repmat(g_x,nrEchoes);
g_x = g_x(1:nrEchoes,1:length(spec2)/nrEchoes); % creates gaussian filter for each echo

gFilt = reshape(g_x',length(spec2),1); %reshapes into a whole echo train 

for i = 1:size(spec2,1) % multiplies the gaussian by the echo train
    spec2(i,:) = spec2(i,:).*gFilt';
end

%% reshaping data
close all
spec3 = reshape(real(spec2'),nrPts, nrEchoes, 360);% reshapes into number of points, number of echoes, total number of exps
spec3sum = sum(spec3);% sums each individual echo for using
spec3sum = reshape (spec3sum, nrEchoes, 8, 45); %echo trains for each of the 45 experiments, 

spec3normSum = spec3sum./max(max(max(spec3sum)));
% 
% for i = 1:45 
%     
tEcho = (1:nrEchoes)*150e-6; % creates time axis
tEcho = tEcho';

spec4 = (spec3normSum(:,:,10)');
save('intdata.txt','spec4','-ascii'); % saves in a text file for use with 2dLaplace
save('tEcho.txt','tEcho','-ascii');% saves in a text file for use with 2dLaplace

% these tD files are taken from the excel sheet

tD1 = [8.6000000e+01
   1.3000000e+02
   1.7400000e+02
   2.1800000e+02
   2.6200000e+02
   3.0600000e+02
   3.5000000e+02
   3.9400000e+02];

tD2 = [8.1000000e+02
   7.2200000e+02
   6.3400000e+02
   5.4600000e+02
   4.5800000e+02
   3.7000000e+02
   2.8200000e+02
   1.9400000e+02];

tD3 = [1.6400000e+02
   2.0800000e+02
   2.5200000e+02
   2.9600000e+02
   3.4000000e+02
   3.8400000e+02
   4.2800000e+02
   4.7200000e+02];

tD1 = tD1./1000;
tD2 = (tD2./1000).^2;
tD3 = tD3./1000;

save('tD1.txt','tD1','-ascii');
save('tD2.txt','tD2','-ascii');
save('tD3.txt','tD3','-ascii');

TwoDLaplaceInverse