clear
close all
clc
addpath(genpath('Z:\TKM\'));

%%Mortar Drying Matlab
maindir = 'C:\Users\bmfortman\Documents\Data\';

dir11 = 'C:\Users\bmfortman\Documents\Data\MortarDrying\OneToOne\1\';
dir21 = 'C:\Users\bmfortman\Documents\Data\MortarDrying\TwoToOne\1\';
dir21T = 'C:\Users\bmfortman\Documents\Data\MortarDrying\TwoToOneThin\2\';
dir21C = 'C:\Users\bmfortman\Documents\Data\MortarDrying\TwoToOneCured\1\';
sample7= 'C:\Users\bmfortman\Documents\Data\MortarDrying\2-1Samp7\1';


%% reading mortar data in then plotting
cd(sample7);

% params.echoTime = readpar_Kea('acqu.par','echoTime');
% params.pulseLength = readpar_Kea('acqu.par','pulseLength');
% params.dwellTime = readpar_Kea('acqu.par','dwellTime');
params.nrExp = readpar_Kea('acqu.par','nrExp');
params.nrSteps = readpar_Kea('acqu.par','nrSteps');
params.nrScans = readpar_Kea('acqu.par','nrScans');
params.repTime = readpar_Kea('acqu.par','repTime');

expLength = params.nrSteps * params.nrExp * params.repTime * params.nrScans/3600000; %for hours /3600000; for minutes /60000; for seconds /1000
time = linspace(0,expLength,(params.nrSteps * params.nrExp)); % this creates a time vector of the length of the data concatenations

expNums = [1 : (params.nrExp)];% don't need to change for the number of experiments

for j = 1:9% this is for experiments with more profiles than 10         
    data = strcat('2-1Samp70',num2str(j),'.dat');% need to change to name of one of the directory file
      

    data1(j).sig = load(data);
end

for j = 10:length(expNums)
        
     data = strcat('2-1Samp7',num2str(j),'.dat'); % need to change to name of one of the directory file for different experiments

    data1(j).sig = load(data);
end

dataAll= [data1(1).sig];

for j = 2:length(expNums)
dataAll = [dataAll; data1(j).sig]; %data1(2).sig; data1(3).sig; data1(4).sig; data1(5).sig; data1(6).sig; data1(7).sig; data1(8).sig;];
end
% expTime = 6.4;
% totTime = expTime*length(dataAll);
% time = linspace(0,totTime,216);% works for dir11, dir21
figure(1) % this gives a 3-d plot
hold on
% plot3(time,(dataAll(:,1)./1000),dataAll(:,2)) % for multiple positions
xlabel('time (minutes)')
ylabel('depth (Millimeters)')
% zlabel('amplitude') % for multiple positions

% expNums = [1 : 9];
% cd(dir21);
% for j = 1:length(expNums)
%     cd(dir21);
%     
%     data = strcat('TwoToOne0',num2str(j),'.dat');
%     data1(j).sig = load(data);
% end
% dataAll2 = [data1(1).sig; data1(2).sig; data1(3).sig; data1(4).sig; data1(5).sig; data1(6).sig; data1(7).sig; data1(8).sig; data1(9).sig];
% % expTime = 6.4;
% % totTime = expTime*length(dataAll2);
% % time = linspace(0,totTime,261);
% figure(1)
% hold on
% plot3(time,dataAll2(:,1),dataAll2(:,2),'-r')
% xlabel('time')
% ylabel('depth')
% zlabel('amplitude')

%% plotting for each depth independently
depth1 = dataAll(1,:); % for each depth you need another depth point
depth2 = dataAll(2,:); 
depth3 = dataAll(3,:);
time1 = time(1);
time2 = time(2);
time3 = time(3);
for j = 4 : (params.nrExp * params.nrSteps) %starting at point 4 since the first 3 are delineated above
    if dataAll(j,1) == 6000 % these are the depths at which the points are measuring at
        depth1 = vertcat(depth1,dataAll(j,:));
        time1 = vertcat(time1,time(j));
    elseif dataAll(j,1) == 3000 % these will need to be changed for different depths
        depth2 = vertcat(depth2,dataAll(j,:));
        time2 = vertcat(time2,time(j));
    else
        depth3 = vertcat(depth3,dataAll(j,:));
        time3 = vertcat(time3,time(j));
    end
  
end

%% gaussian fitting
guesses = [1.2, 7000];
%f = fit(time2,(depth2(:,2)),'gauss1');
[fit(1).beta,fit(1).resid,fit(1).J] = nlinfit(time1,(depth1(:,2)),@gauss_simple,guesses);
[fit(2).beta,fit(2).resid,fit(2).J] = nlinfit(time2,(depth2(:,2)),@gauss_simple,guesses);
[fit(3).beta,fit(3).resid,fit(3).J] = nlinfit(time3,(depth3(:,2)),@gauss_simple,guesses);

fit(1).pred = gauss_simple(fit(1).beta,(time1));
fit(2).pred = gauss_simple(fit(2).beta,(time2));
fit(3).pred = gauss_simple(fit(3).beta,(time3));

figure(2)
hold on
plot(time1,(depth1(:,2)))
plot(time1,(fit(1).pred))
plot(time2,(depth2(:,2)),'-r')
plot(time2,(fit(2).pred),'-r')
plot(time3,(depth3(:,2)),'-k')
plot(time3,(fit(3).pred),'-k')
xlabel('Time (Seconds)') %this depends upon what you selected the time to be displayed in line 26
ylabel('Amplitude')
legend('6mm depth (bottom)','3mm depth','0mm depth (top)')% edit these if needed, with appropriate depths from above

% figure(3)
% hold on
% plot(time1.^2,(depth1(:,2)))
% plot(time2.^2,(depth2(:,2)),'-r')
% plot(time3.^2,(depth3(:,2)),'-k')
% xlabel('Time Squared (Seconds^2)') %this depends upon what you selected the time to be displayed in line 26
% ylabel('Amplitude')
% legend('6mm depth (bottom)','3mm depth','0mm depth (top)')% edit these if needed, with appropriate depths from above



%% reading t&H data in then plotting
dirTH = 'C:\Users\bmfortman\Documents\Data\T&HDatalogger\';
cd(maindir)
cd(dirTH);

OneToOneTH = csvread('MortarDryingOneToOneEdit.txt',0,2);%this is in C
TwoToOneTH = csvread('MortarDryingTwoToOneEdit.txt',0,2); %this is in F needs to be converted

Temp21c = ((TwoToOneTH(:,1) - 32)/1.8);

figure(4)
hold on
axis([630 1980 20 35]);
plot(OneToOneTH(:,1))%this is temp
plot(OneToOneTH(:,2),'-r')%this is humidityin %rh
legend('Temp1:1 Deg.C ','Humidity1:1 % humidity')

figure(5)
hold on
axis([28 1010 20 30])
plot(Temp21c,'-g')
plot(TwoToOneTH(:,2),'-k');
legend('Temp2:1 Deg C','Humidity2:1 % humidity')