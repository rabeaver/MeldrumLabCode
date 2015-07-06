close all
clear
clc
%% input Data file
cd('C:\Users\benjamin\Documents\Data\Mortar Curing\CO2 meter test experiments')

data = csvread('November-04-20142.csv',0,1);
data = (data(316:end,1));
data = data'./max(data);
timeVector = ((1:length(data))*0.5)/60;% in minutes

% % exponential fit
guesses = [0.9,-0.003];
[beta,r,j] = nlinfit(timeVector,data,@exponentialFxn,guesses);
Sfit = exponentialFxn(beta,timeVector);

% % linear fit
% [P] = polyfit(timeVector,data,1);
% Sfit = P(1)*timeVector +P(2);

%% plot
figure(1)
plot(timeVector,data);
hold on
plot(timeVector,Sfit);
