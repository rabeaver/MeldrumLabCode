% close all
clear
clc
%% setup of params and dir
labCompy = ('C:\Users\bmfortman\Documents\Data\');
sample1_1 = ('MortarDrying\OneToOne\1\');
sample2_1 = ('MortarDrying\TwoToOne\1');

cd(strcat(labCompy,sample2_1));

params.nrExp = readpar_Kea('acqu.par','nrExp');
params.nrSteps = readpar_Kea('acqu.par','nrSteps');
params.nrScans = readpar_Kea('acqu.par','nrScans');
params.repTime = readpar_Kea('acqu.par','repTime');

expLength = params.nrSteps * params.nrExp * params.repTime * params.nrScans/3600000; %for hours /3600000; for minutes /60000; for seconds /1000
time = linspace(0,expLength,11); % this creates a time vector of the length of the data concatenations

expNums = [1 : (params.nrExp)];% don't need to change for the number of experiments

%% read in data
for i = 1%:nrExp;
    data = load(strcat('TwoToOne0',num2str(i),'.dat'));
end

for i = 2:params.nrExp;
    data2 = load(strcat('TwoToOne0',num2str(i),'.dat'));
    data = [data,data2(:,2)];
end

depthAxis = data(:,1)./10000;
data = data(:,2:end);
data = data./max(max(data));

%% plotting data
figure(1)
hold on
surf(time,depthAxis,data)
shading flat
xlabel('Time [hours]')
ylabel('Depth [cm]')
zlabel('Relative Amplitude')
title('One to One Drying')

for j = 1:size(data,2)
timeN(j).pts = ones(1,size(data,1)).*time(j);
end

timeNew = timeN(1).pts;
depthNew = depthAxis;
dataNew = reshape(data,1,(size(data,1)*size(data,2)));
for i = 2:size(data,2)
    timeNew = [timeNew,timeN(i).pts];
    depthNew = [depthNew;depthAxis];
end
figure(2)
hold on
plot3(timeNew,depthNew,dataNew,'k','LineWidth',2)
xlabel('Time [hours]')
ylabel('Depth [cm]')
zlabel('Relative Amplitude')