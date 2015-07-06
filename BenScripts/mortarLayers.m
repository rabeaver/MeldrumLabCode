clear
close all
clc
addpath(genpath('Z:\TKM\'));

%% Mortar layers Matlab
maindir = 'C:\Users\bmfortman\Documents\Data\';

dir1 = 'C:\Users\bmfortman\Documents\Data\Layers\Anchoring\3';
dir2 = 'C:\Users\bmfortman\Documents\Data\Layers\QuickSet\1';
cd(maindir)

%% loading data
cd(dir1)
aData = load('Anchoring.dat');

cd(dir2)
qData = load('QuickSet.dat');

%% smoothing data
% aSmooth1 = sgolayfilt(aData(:,2),3,11);
% aSmooth2 = sgolayfilt(aData(:,2),4,11);
% aSmooth3 = sgolayfilt(aData(:,2),5,11);
% aSmooth4 = sgolayfilt(aData(:,2),6,11);
% arSmooth1 = sgolayfilt(aData(:,2),5,7);
% arSmooth2 = sgolayfilt(aData(:,2),5,11);
% arSmooth3 = sgolayfilt(aData(:,2),5,13);
arSmooth4 = sgolayfilt(aData(:,2),5,15);
qSmooth = sgolayfilt(qData(:,2),5,15);
derA = diff(arSmooth4)./diff(aData(:,1));
derQ = diff(qSmooth)./diff(qData(:,1));

figure(1)
hold on
plot(aData(:,1),aData(:,2))
%plot(aData(:,1),aSmooth1,'-r')
% plot(aData(:,1),aSmooth2,'-g')
% plot(aData(:,1),aSmooth3,'-k')
% plot(aData(:,1),aSmooth4,'-y')
% 
% figure(2)
% hold on
% plot(aData(:,1),aData(:,2))
%plot(aData(:,1),arSmooth1,'-y')
%plot(aData(:,1),arSmooth2,'-g')
%plot(aData(:,1),arSmooth3,'-k')
plot(aData(:,1),arSmooth4,'-k')
legend('Anchoring Data', 'Anchoring Smoothed')
xlabel('postion')
ylabel('Signal Amplitude')

figure(2) 
hold on
plot(qData(:,1),qData(:,2))
plot(qData(:,1),qSmooth,'-k')
legend('Quickset Data','Quickset Smoothed')
xlabel('position')
ylabel('Signal Amplitude')

figure(3)
hold on
axis ([0 14000 -1e-4 1e-4])
plot(aData((2:end),1),derA)
plot(qData((2:end),1),derQ,'-r')
legend('Derivative of Anchoring','Derivative of Quickset')
xlabel('position')
