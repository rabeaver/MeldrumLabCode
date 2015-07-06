close all
clear all
clc

cd('C:\Users\benjamin\Documents\Data\MortarDrying\1-1Samp1\1');
% these are all for sample 40 on my laptop, saved using the commented out
% section on lines 167-179 of the
% mortardryingt2ampsinglelayersampleanalysis script
%% loading data from the saved txt files
%name of the samples
savename(1) = ('Armory')
savename(2) = ('Jamestown')
savename(3) = ('Quickset')
savename(4) =

dataNorm = load('dataNorm.txt');
echoAxis = load('echoAxis.txt');

dataPredBit2fix = load('dataPredBit2fix.txt');
dataResidBit2fix = load('dataResidBit2fix.txt');

dataPredBi = load('dataPredBi.txt');
dataResidBi = load('dataResidBi.txt');

dataPredMono = load('dataPredMono.txt');
dataResidMono = load ('dataResidMono.txt');

zeroLine = zeros(length(echoAxis),1); %creates a line of zeroes for comparison with residuals
%% plotting the data into one figure
figure(1)% normalized data versus the fits
hold on 
plot(echoAxis,dataNorm,'-m');
plot(echoAxis,dataPredBit2fix,'-k');
plot(echoAxis,dataPredBi,'-r');
plot(echoAxis,dataPredMono,'-b');
xlabel('time [us]')
ylabel('Normalized summed Echoes')
legend('Raw data','Predicted Biexp fit Fixed 2nd T2','Predicted Biexp fit','Predicted MonoExp fit');

figure(2)
hold on
plot(echoAxis,dataResidBit2fix,'-k');
axis([0 max(echoAxis) -0.03 0.05])
plot(echoAxis,dataResidBi,'-r');
plot(echoAxis,dataResidMono,'-b');
plot(echoAxis,zeroLine,'-k','LineWidth',2);
xlabel('time [us]')
title('Residuals from fitting')
legend('Predicted Biexp fit Fixed 2nd T2','Predicted Biexp fit','Predicted MonoExp fit','Zeroline');

