clear
close all
clc

cd '/Users/tyler/Dropbox/Data/Biosensors/';
load('T2D.mat');

figure(1)
subplot(1,2,2)
hold on
plot(T2(:,1),D(:,1),'-ok')
plot(T2(:,2),D(:,2),'-or')
text(T2(1,1)*1.1, D(1,1), 'Point 1');
xlim([0 100])
ylim([0 1.8e-10])
xlabel('T_2 [ms]')
ylabel('D [m^2 s^{-1}]')
subplot(1,2,1)
hold on
plot(T2(:,1),D(:,1),'-ok')
plot(T2(:,2),D(:,2),'-or')
text(T2(2,2)*1.3, D(2,2), 'Point 2');
xlim([0 5])
ylim([0 1.8e-10])
xlabel('T_2 [ms]')
ylabel('D [m^2 s^{-1}]')