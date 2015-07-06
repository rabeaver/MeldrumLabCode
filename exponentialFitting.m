clear
clc
close all

dirhead = 'C:\Users\TMeldrum\Dropbox\Data\Vinh\Paint Data\5.10.12\crimson\'; % data dir
datafile = 'cpmg5_10_12-500-512.txt'; % data file name
CI = 90; %confidence interval (in percent, ie '90' = 90% confidence
data = load(strcat(dirhead,datafile));

%% Monoexponential fit
y0_guess = 0; %initial guess, y0
A_guess = 10; %initial guess, A
t2_guess = 14.6; %initial guess, tau


guesses = [y0_guess;A_guess;t2_guess]; 
[xfitmono,ypredmono,coeffsmono,residualsmono] = monodecay_t2fit(data(:,1),data(:,2),guesses,CI);

%% Biexponential fit
y0_guess = 0; %initial guess, y0
A_guess1 = 10; %initial guess, A1
t2_guess1 = 33.4; %initial guess, tau1
A_guess2 = 5; %initial guess, A2
t2_guess2 = 2.1; %initial guess, tau2

guesses = [y0_guess;A_guess1;t2_guess1;A_guess2;t2_guess2];
[xfitbi,ypredbi,coeffsbi,residualsbi] = bidecay_t2fit(data(:,1),data(:,2),guesses,CI);

%% Plotting
figure
subplot(2,1,1)
hold on
scatter(data(:,1),data(:,2),'ob')
plot(xfitmono,ypredmono,'-k')
plot(xfitbi,ypredbi,'-r')
legend('data','monoexponential fit','biexponential fit')
subplot(2,1,2)
hold on
bar(data(:,1),[residualsmono,residualsbi])
legend('monoexp. residuals','biexp. residuals')