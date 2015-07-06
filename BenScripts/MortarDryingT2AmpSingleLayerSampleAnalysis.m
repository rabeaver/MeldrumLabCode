clear
close all
clc
% addpath(genpath('Z:\TKM\'));


%% Mortar Drying Matlab
compStr = 'C:\Users\benjamin\Documents\Data'; % for my laptop
% compStr = 'C:\Users\bmfortman\Documents\Data';%for the lab compy

sample(1).dir = strcat(compStr,'\MortarDrying\1-1Samp1\1');% combines the computer string with the name of the sample so that it runs for 
sample(2).dir = strcat(compStr,'\MortarDrying\Sample2\1');
sample(3).dir = strcat(compStr,'\MortarDrying\Sample3\1');
sample(4).dir = strcat(compStr,'\MortarDrying\Sample4\1');
sample(5).dir = strcat(compStr,'\MortarDrying\Sample5\1');
sample(6).dir = strcat(compStr,'\MortarDrying\Sample6\1');
sample(7).dir = strcat(compStr,'\MortarDrying\2-1Samp7\1');
sample(8).dir = strcat(compStr,'\MortarDrying\Sample8\1');
sample(9).dir = strcat(compStr,'\MortarDrying\Sample9\1');
sample(10).dir = strcat(compStr,'\MortarDrying\Sample10\1');
sample(11).dir = strcat(compStr,'\MortarDrying\2-1Samp11\1');
sample(12).dir = strcat(compStr,'\MortarDrying\Sample12\1');

sample(92).dir = strcat(compStr,'\MortarDrying\Sample92\1'); % these are the retesting samples 2-1
sample(112).dir = strcat(compStr,'\MortarDrying\Sample112\1');% 2-1
sample(13).dir = strcat(compStr,'\MortarDrying\Sample13\1'); % 1-1 1% brick dust
sample(52).dir = strcat(compStr,'\MortarDrying\Sample52\1');% this IS actually a rerun of sample 5, the only difference being time before running
sample(72).dir = strcat(compStr,'\MortarDrying\Sample72\1');

sampleNums = [1 2 3 4 5 6 7 8 9 10 11 12 13 52 92 112 72]; %these are the sample numbers that have been run, will need to be updated with new ones and a matching directory, and a plot color in the final figure
for i = sampleNums; 
cd(sample(i).dir); %when changing directories you also need to change the name of the file properly in lines 30,31,38,39

%% Reading in Paramaters

params.echoTime = readpar_Kea('acqu.par','echoTime');
params.pulseLength = readpar_Kea('acqu.par','pulseLength');
params.dwellTime = readpar_Kea('acqu.par','dwellTime');
params.nrExp = readpar_Kea('acqu.par','nrExp');
params.nrSteps = readpar_Kea('acqu.par','nrSteps');
params.nrScans = readpar_Kea('acqu.par','nrScans');
params.repTime = readpar_Kea('acqu.par','repTime');

%% reading mortar data in then plotting
expNums = [1 : (params.nrExp)]; % should cover the number of experiments, if the experiment was cut short, simply change the #exp in the acqu.par file before running
if i == 1 % the purpose of this loop is to change the reading parameters for the experiments that have different names than sample 5
    
for j = 1:9% works provided there are more than 9 experiments
        
    dataR = strcat('2-1Samp10',num2str(j),'-decaysRe.dat');
    dataI = strcat('2-1Samp10',num2str(j),'-decaysIm.dat');
    dataA = strcat('2-1Samp10',num2str(j),'.dat');
    data(j).amp = load(dataA);
    data(j).real = load(dataR);
    data(j).imag = load(dataI);
end

for j = 10:length(expNums) % covers the rest of the experiments
        
    dataR = strcat('2-1Samp1',num2str(j),'-decaysRe.dat'); % need to change to name of one of the directory file for different experiments
    dataI = strcat('2-1Samp1',num2str(j),'-decaysIm.dat');
    dataA = strcat('2-1Samp1',num2str(j),'.dat');
    data(j).amp = load(dataA);
    data(j).real = load(dataR);
    data(j).imag = load(dataI);
end

elseif i == 7 %covers experiment 7
    
for j = 1:9% works provided there are more than 9 experiments
        
    dataR = strcat('2-1Samp70',num2str(j),'-decaysRe.dat');
    dataI = strcat('2-1Samp70',num2str(j),'-decaysIm.dat');
    dataA = strcat('2-1Samp70',num2str(j),'.dat');
    data(j).amp = load(dataA);
    data(j).real = load(dataR);
    data(j).imag = load(dataI);
end

for j = 10:length(expNums) % covers the rest of the experiments
        
    dataR = strcat('2-1Samp7',num2str(j),'-decaysRe.dat'); % need to change to name of one of the directory file for different experiments
    dataI = strcat('2-1Samp7',num2str(j),'-decaysIm.dat');
    dataA = strcat('2-1Samp7',num2str(j),'.dat');
    data(j).amp = load(dataA);
    data(j).real = load(dataR);
    data(j).imag = load(dataI);
end
elseif i == 11 %covers experiment 11
    
for j = 1:9% works provided there are more than 9 experiments
        
    dataR = strcat('2-1Samp110',num2str(j),'-decaysRe.dat');
    dataI = strcat('2-1Samp110',num2str(j),'-decaysIm.dat');
    dataA = strcat('2-1Samp110',num2str(j),'.dat');
    data(j).amp = load(dataA);
    data(j).real = load(dataR);
    data(j).imag = load(dataI);
end

for j = 10:length(expNums) % covers the rest of the experiments
        
    dataR = strcat('2-1Samp11',num2str(j),'-decaysRe.dat'); % need to change to name of one of the directory file for different experiments
    dataI = strcat('2-1Samp11',num2str(j),'-decaysIm.dat');
    dataA = strcat('2-1Samp11',num2str(j),'.dat');
    data(j).amp = load(dataA);
    data(j).real = load(dataR);
    data(j).imag = load(dataI);
end
    
else %this will cover all experiments with the name 'Sample(i)' where i is the number of the sample entered in the directory above
    
for j = 1:9% works provided there are more than 9 experiments
        
    dataR = strcat('Sample',num2str(i),'0',num2str(j),'-decaysRe.dat');
    dataI = strcat('Sample',num2str(i),'0',num2str(j),'-decaysIm.dat');
    dataA = strcat('Sample',num2str(i),'0',num2str(j),'.dat');
    data(j).amp = load(dataA);
    data(j).real = load(dataR);
    data(j).imag = load(dataI);
end

for j = 10:length(expNums) % covers the rest of the experiments
        
    dataR = strcat('Sample',num2str(i),num2str(j),'-decaysRe.dat'); % need to change to name of one of the directory file for different experiments
    dataI = strcat('Sample',num2str(i),num2str(j),'-decaysIm.dat');
    dataA = strcat('Sample',num2str(i),num2str(j),'.dat');
    data(j).amp = load(dataA);
    data(j).real = load(dataR);
    data(j).imag = load(dataI);
end
    
end

echoAxis = [data(1).real(3:end,1)]; % works for the kea, not so for the tecMag
dataReal = [data(1).real(3:end,2), data(1).real(3:end,3)]; % these are starting at the 3rd point to avoid extra noise from the pulse sequence
dataImag = [data(1).imag(3:end,2), data(1).imag(3:end,3)];% the 2nd and 3rd columns are depths, for samples here it is all at the same depth so these are disregarded
dataAmp = [data(1).amp(:,2)]; % takes the amplitudes for the experiment

for j = 2:length(expNums) % THIS puts all of the data into a single matrix 
    dataReal = [dataReal, data(j).real(3:end,2), data(j).real(3:end,3)];
    dataImag = [dataImag, data(j).imag(3:end,2), data(j).imag(3:end,3)];
    dataAmp = [dataAmp;data(j).amp(:,2)];
    
end
%% normalizing Data
maxS = max(dataReal);
maxA = max(dataAmp);
dataNorm = dataReal(:,1)/maxS(1);
sample(i).ampNorm = dataAmp(:,1)/maxA(1);
for j = 2:size(dataReal,2);
    dataNorm = [dataNorm, dataReal(:,j)/maxS(j)];
end
%% exponential fit

guesses = [0.9,18,0.1];% [1, 15];% guesses for the nlin fit
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');


for j = 1:size(dataNorm,2);% this is to move along the matrices taking the next point
    fit(j).beta=[0,0,0];
    while fit(j).beta~=guesses %generates while loop if guess is off
      

    [fit(j).beta,fit(j).resid,fit(j).J] = nlinfit(echoAxis,dataNorm(:,j),@t2bifit_t22fixed,guesses,fitopts);% exp fit bifit__t22fixed has a fixed value of 4.42us
    fit(j).pred = t2bifit_t22fixed(fit(j).beta,echoAxis); % produces predicted values for easy visualization of fit
%     for j = 40;
%     figure(j) % used for visualizing the fit and saving data from fits
%     hold on % if necessary
%     plot(echoAxis,dataNorm(:,j))
%     plot(echoAxis,fit(j).pred,'-k')
%     plot(echoAxis,fit(j).resid,'-r')
%     Norm = dataNorm(:,j);
%     Pred = fit(j).pred;
%     Resid = fit(j).resid;
%     save('dataNorm.txt','Norm','-ascii')
%     save('dataPredBit2fix.txt','Pred','-ascii')
%     save('dataResidBit2fix.txt','Resid','-ascii')
%     end
    guesses=fit(j).beta; %reallocates guess to new result
    end
    ci = nlparci(fit(j).beta,fit(j).resid,'jacobian',fit(j).J);
    error(j).ci = ci; % gives confidence intervals
  
end

for j = 1:size(dataNorm,2); % converts the structures to matrices for ease of graphing
    exp(j) = fit(j).beta(1);
    t2(j) = fit(j).beta(2);
    eExp(1,j) = error(j).ci(1);
    eExp(2,j) = error(j).ci(3);
    eT2(1,j) = error(j).ci(2);
    eT2(2,j) = error(j).ci(4);
    
    exp2(j) = fit(j).beta(3); %These are for the biexp fit
    % t22(j) = fit(j).beta(4);
    eExp2(1,j) = error(j).ci(5);
    eExp2(2,j) = error(j).ci(6);
%     eT22(1,j) = error(j).ci(7);
%     eT22(2,j) = error(j).ci(8);
    
    
%     resid(j) = fit(j).resid(:);
%     pred(j) = fit(j).pred(:);
end

sample(i).exp = exp;
sample(i).t2 = t2;
sample(i).eExpdiffl =  abs(eExp(1,:)-exp); %converts confidence intervals into distance from values for use with the errorbar function
sample(i).eExpdiffu =  (eExp(2,:)-exp);
sample(i).eT2diffl =  abs(eT2(1,:)-t2);
sample(i).eT2diffu =  (eT2(1,:)-t2);

sample(i).exp2 = exp2;
% sample(i).t22 = t22;
% sample(i).t22Avg = sum(t22)/size(sample(i).t22,2);
% 
% sample(i).resid = fit.resid;
% sample(i).pred = fit.pred;
% sample(i).dataNorm = dataNorm;

%% plotting Monoexpfit vs. time and depth
close all

expTime = params.repTime * params.nrScans * params.nrExp * params.nrSteps/3600000;% for seconds/1000; for minutes/60000; for hours /3600000
sample(i).timePoints = linspace(0,expTime,size(exp,2));
clear('maxS','dataReal','dataImag','data','dataNorm','eExp','eT2','error','echoAxis','t2','fit','exp','resid','dataAmp','exp2');
end

figure(1)
subplot(2,1,1);
hold on
plot(sample(1).timePoints,sample(1).exp,'-b')% these will just plot the exponents
axis([0 28 -0.1 1.2])
plot(sample(2).timePoints,sample(2).exp,'-m')
plot(sample(3).timePoints,sample(3).exp,'-c')
plot(sample(4).timePoints,sample(4).exp,'-r')
plot(sample(5).timePoints,sample(5).exp,'-k')
plot(sample(6).timePoints,sample(6).exp,'-g')
plot(sample(13).timePoints,sample(13).exp,'-c','LineStyle','--','LineWidth',1.5)

plot(sample(1).timePoints,sample(1).exp2,'-b')% these will just plot the exponents
plot(sample(2).timePoints,sample(2).exp2,'-m')
plot(sample(3).timePoints,sample(3).exp2,'-c')
plot(sample(4).timePoints,sample(4).exp2,'-r')
plot(sample(5).timePoints,sample(5).exp2,'-k')
plot(sample(6).timePoints,sample(6).exp2,'-g')
plot(sample(13).timePoints,sample(13).exp2,'-c','LineStyle','--','LineWidth',1.5)


% errorbar(sample(1).timePoints,sample(1).exp,sample(1).eExpdiffl,sample(1).eExpdiffu)% errorbars are a bit messy, so giving option to remove them
% errorbar(sample(5).timePoints,sample(5).exp,sample(5).eExpdiffl,sample(5).eExpdiffu,'-r')
% errorbar(sample(72).timePoints,sample(72).exp,sample(72).eExpdiffl,sample(72).eExpdiffu,'-k')
% errorbar(sample(11).timePoints,sample(11).exp,sample(11).eExpdiffl,sample(11).eExpdiffu,'-g')
% errorbar(sample(2).timePoints,sample(2).exp,sample(2).eExpdiffl,sample(2).eExpdiffu,'-m')
title('One to one lime to sand samples')
xlabel('Experiment Time [hours]')
ylabel('Exponent')
legend('Sample1 No Additives','Sample2 Ash Only','Sample3 Brick Dust Only','Sample4 Clay Only','Sample5 All Additives','Sample6 Clay/Brick Dust','Sample13 1% brick dust')%, '2nd Biexponent')

subplot(2,1,2);
hold on
plot(sample(72).timePoints,sample(72).exp,'-b')% these will just plot the Amplitudes
axis([0 28 -0.1 1.2])
plot(sample(8).timePoints,sample(8).exp,'-m')
plot(sample(92).timePoints,sample(92).exp,'-c')
plot(sample(10).timePoints,sample(10).exp,'-r')
plot(sample(112).timePoints,sample(112).exp,'-k')
plot(sample(12).timePoints,sample(12).exp,'-g')

plot(sample(72).timePoints,sample(72).exp2,'-b')% these will just plot the Amplitudes
plot(sample(8).timePoints,sample(8).exp2,'-m')
plot(sample(92).timePoints,sample(92).exp2,'-c')
plot(sample(10).timePoints,sample(10).exp2,'-r')
plot(sample(112).timePoints,sample(112).exp2,'-k')
plot(sample(12).timePoints,sample(12).exp2,'-g')


% errorbar(sample(1).timePoints,sample(1).exp,sample(1).eExpdiffl,sample(1).eExpdiffu)% errorbars are a bit messy, so giving option to remove them
% errorbar(sample(5).timePoints,sample(5).exp,sample(5).eExpdiffl,sample(5).eExpdiffu,'-r')
% errorbar(sample(72).timePoints,sample(72).exp,sample(72).eExpdiffl,sample(72).eExpdiffu,'-k')
% errorbar(sample(11).timePoints,sample(11).exp,sample(11).eExpdiffl,sample(11).eExpdiffu,'-g')
% errorbar(sample(2).timePoints,sample(2).exp,sample(2).eExpdiffl,sample(2).eExpdiffu,'-m')
title('Two to one lime to sand samples')
xlabel('Experiment Time [hours]')
ylabel('Exponent')
legend('Sample7 No Additives','Sample8 Ash Only','Sample9 Brick Dust Only','Sample10 Clay Only','Sample11 All Additives','Sample12 Clay/Brick Dust')%, '2nd Biexponent')

figure(2) % this will plot the normalized amplitudes so that the drying curve can easily be seen
subplot(1,2,1);
hold on
axis([0 25 0 1.1])
plot(sample(1).timePoints,sample(1).ampNorm)% these will just plot the exponents
plot(sample(2).timePoints,sample(2).ampNorm,'-m')
plot(sample(3).timePoints,sample(3).ampNorm,'-c')
plot(sample(4).timePoints,sample(4).ampNorm,'-r')
plot(sample(5).timePoints,sample(5).ampNorm,'-k')
plot(sample(6).timePoints,sample(6).ampNorm,'-g')
plot(sample(13).timePoints,sample(13).ampNorm,'-c','LineWidth',1.5,'LineStyle','--')

title('50% lime 50% sand samples')
xlabel('Experiment Time [hours]')
ylabel('Normalized Amplitude')
legend('Sample1 No Additives','Sample2 Ash Only','Sample3 Brick Dust Only','Sample4 Clay Only','Sample5 All Additives','Sample6 Clay/Brick Dust','Sample13 1% brick Dust')%, '2nd Biexponent')

% figure(15)
subplot(1,2,2);
hold on
axis([0 25 0 1.1])
plot(sample(72).timePoints,sample(72).ampNorm)% these will just plot the Amplitudes
plot(sample(8).timePoints,sample(8).ampNorm,'-m')
plot(sample(92).timePoints,sample(92).ampNorm,'-c')
plot(sample(10).timePoints,sample(10).ampNorm,'-r')
plot(sample(112).timePoints,sample(112).ampNorm,'-k')
plot(sample(12).timePoints,sample(12).ampNorm,'-g')


% errorbar(sample(1).timePoints,sample(1).exp,sample(1).eExpdiffl,sample(1).eExpdiffu)% errorbars are a bit messy, so giving option to remove them
% errorbar(sample(5).timePoints,sample(5).exp,sample(5).eExpdiffl,sample(5).eExpdiffu,'-r')
% errorbar(sample(72).timePoints,sample(72).exp,sample(72).eExpdiffl,sample(72).eExpdiffu,'-k')
% errorbar(sample(11).timePoints,sample(11).exp,sample(11).eExpdiffl,sample(11).eExpdiffu,'-g')
% errorbar(sample(2).timePoints,sample(2).exp,sample(2).eExpdiffl,sample(2).eExpdiffu,'-m')
title('66% lime 33% sand samples')
xlabel('Experiment Time [hours]')
ylabel('Normalized Amplitude')
legend('Sample7 No Additives','Sample8 Ash Only','Sample9 Brick Dust Only','Sample10 Clay Only','Sample11 All Additives','Sample12 Clay/Brick Dust')%, '2nd Biexponent')



figure(3)
hold on
axis([0 20 5 30])
errorbar(sample(1).timePoints,sample(1).t2,sample(1).eT2diffl,sample(1).eT2diffu)
errorbar(sample(2).timePoints,sample(2).t2,sample(2).eT2diffl,sample(2).eT2diffu,'-m')
errorbar(sample(3).timePoints,sample(3).t2,sample(3).eT2diffl,sample(3).eT2diffu,'-c')
errorbar(sample(4).timePoints,sample(4).t2,sample(4).eT2diffl,sample(4).eT2diffu,'-r')
errorbar(sample(5).timePoints,sample(5).t2,sample(5).eT2diffl,sample(5).eT2diffu,'-k')
errorbar(sample(6).timePoints,sample(6).t2,sample(6).eT2diffl,sample(6).eT2diffu,'-g')
errorbar(sample(13).timePoints,sample(13).t2,sample(13).eT2diffl,sample(13).eT2diffu,'-c','LineWidth',1.5,'LineStyle','--')

% plot(sample(1).timePoints,-1.114 * (sample(1).timePoints) + 25.78)%
% plot(sample(2).timePoints,-1.219 * (sample(2).timePoints) + 26.41,'-m')
% plot(sample(3).timePoints,-1.159 * (sample(3).timePoints) + 24.64,'-c')
% plot(sample(4).timePoints,-1.37 * (sample(4).timePoints) + 24.97,'-r')
% plot(sample(5).timePoints,-1.112 * (sample(5).timePoints) + 23.49,'-k')
% plot(sample(6).timePoints,-1.292 * (sample(6).timePoints) + 23.25,'-g')

%plot(expDepth,biExpt22,'-r')
title('One to One Samples')
xlabel('Experiment Time [hours]')
ylabel('T2 time [ms]')
legend('Sample1 No Additives','Sample2 Ash Only','Sample3 Brick Dust Only','Sample4 Clay Only','Sample5 All Additives','Sample6 Clay/Brick Dust','Sample 13 1% Brick Dust')%, '2nd T2 time')

figure(4)
hold on
axis([0 20 5 30])
errorbar(sample(72).timePoints,sample(72).t2,sample(72).eT2diffl,sample(72).eT2diffu)
errorbar(sample(8).timePoints,sample(8).t2,sample(8).eT2diffl,sample(8).eT2diffu,'-m')
errorbar(sample(92).timePoints,sample(92).t2,sample(92).eT2diffl,sample(92).eT2diffu,'-c')
errorbar(sample(10).timePoints,sample(10).t2,sample(10).eT2diffl,sample(10).eT2diffu,'-r')
errorbar(sample(112).timePoints,sample(112).t2,sample(112).eT2diffl,sample(112).eT2diffu,'-k')
errorbar(sample(12).timePoints,sample(12).t2,sample(12).eT2diffl,sample(12).eT2diffu,'-g')

% plot(sample(72).timePoints,-0.5837 * (sample(72).timePoints) + 21.28)% sample slopes calculated in cftool
% plot(sample(8).timePoints,-1.219 * (sample(8).timePoints) + 28.4,'-m')
% plot(sample(9).timePoints,-1.144 * (sample(9).timePoints) + 24.67,'-c')
% plot(sample(10).timePoints,-0.8942 * (sample(10).timePoints) + 26.89,'-r')
% plot(sample(11).timePoints,-0.6193 * (sample(11).timePoints) + 26.38,'-k')
% plot(sample(12).timePoints,-0.8452 * (sample(12).timePoints) + 22.79,'-g')

%plot(expDepth,biExpt22,'-r')
title('Two to One Samples')
xlabel('Experiment Time (hours)')
ylabel('T2 time [ms]')
legend('Sample7 No Additives','Sample8 Ash Only','Sample9 Brick Dust Only','Sample10 Clay Only','Sample11 All Additives','Sample12 Clay/Brick Dust')%, '2nd T2 time')

%% creating derivatives of the T2 Fxn's

for i = sampleNums
    % interpolating to create smoother derivatives
    sample(i).qp = linspace(0,max(sample(i).timePoints),26); %creating # query points from 0 to the max of the time vector Playing around with the number of points to give a good estimation of what is happening
    sample(i).intp = interp1(sample(i).timePoints(:),sample(i).t2(:),sample(i).qp,'spline'); %Creating the interpolation points for the # of querypoints

    % creating the derivatives
    sample(i).t2Deriv = diff(sample(i).t2)./diff(sample(i).timePoints); % derivatives of the t2's
    sample(i).t2DerivIntp = diff(sample(i).intp)./diff(sample(i).qp); % derivatives of the interpolated points
end

% figure(5)
% hold on
% axis([0 20 -2 1])
% plot(sample(1).timePoints(2:end),sample(1).t2Deriv)
% plot(sample(2).timePoints(2:end),sample(2).t2Deriv,'-m')
% plot(sample(3).timePoints(2:end),sample(3).t2Deriv,'-c')
% plot(sample(4).timePoints(2:end),sample(4).t2Deriv,'-r')
% plot(sample(5).timePoints(2:end),sample(5).t2Deriv,'-k')
% plot(sample(6).timePoints(2:end),sample(6).t2Deriv,'-g')
% 
% title('One to One Samples')
% legend('Sample1 No Additives','Sample2 Ash Only','Sample3 Brick Dust Only','Sample4 Clay Only','Sample5 All Additives','Sample6 Clay/Brick Dust')
% title('T2 Derivatives')
% xlabel('Experiment Time (hours)')
% 
% figure(6)
% hold on
% axis([0 20 -2 1])
% plot(sample(72).timePoints(2:end),sample(72).t2Deriv)
% plot(sample(8).timePoints(2:end),sample(8).t2Deriv,'-m')
% plot(sample(9).timePoints(2:end),sample(9).t2Deriv,'-c')
% plot(sample(10).timePoints(2:end),sample(10).t2Deriv,'-r')
% plot(sample(11).timePoints(2:end),sample(11).t2Deriv,'-k')
% plot(sample(12).timePoints(2:end),sample(12).t2Deriv,'-g')
% 
% title('Two to One Samples')
% legend('Sample7 No Additives','Sample8 Ash Only','Sample9 Brick Dust Only','Sample10 Clay Only','Sample11 All Additives','Sample12 Clay/Brick Dust')
% title('T2 Derivatives')
% xlabel('Experiment Time (hours)')

figure(7)
hold on
axis([0 20 -3 1])
plot(sample(1).qp(2:end),sample(1).t2DerivIntp)
plot(sample(2).qp(2:end),sample(2).t2DerivIntp,'-m')
plot(sample(3).qp(2:end),sample(3).t2DerivIntp,'-c')
plot(sample(4).qp(2:end),sample(4).t2DerivIntp,'-r')
plot(sample(5).qp(2:end),sample(5).t2DerivIntp,'-k')
plot(sample(6).qp(2:end),sample(6).t2DerivIntp,'-g')
plot(sample(13).qp(2:end),sample(13).t2DerivIntp,'-c','LineWidth',1.5,'LineStyle','--')


title('One to One Samples')
legend('Sample1 No Additives','Sample2 Ash Only','Sample3 Brick Dust Only','Sample4 Clay Only','Sample5 All Additives','Sample6 Clay/Brick Dust','Sample 13 1% Brick Dust')
title('Interpolated T2 Derivatives')
xlabel('Experiment Time (hours)')

figure(8)
hold on
axis([0 20 -3 1])
plot(sample(72).qp(2:end),sample(72).t2DerivIntp)
plot(sample(8).qp(2:end),sample(8).t2DerivIntp,'-m')
plot(sample(92).qp(2:end),sample(92).t2DerivIntp,'-c')
plot(sample(10).qp(2:end),sample(10).t2DerivIntp,'-r')
plot(sample(112).qp(2:end),sample(112).t2DerivIntp,'-k')
plot(sample(12).qp(2:end),sample(12).t2DerivIntp,'-g')

title('Two to One Samples')
legend('Sample7 No Additives','Sample8 Ash Only','Sample9 Brick Dust Only','Sample10 Clay Only','Sample11 All Additives','Sample12 Clay/Brick Dust')
title('Interpolated T2 Derivatives')
xlabel('Experiment Time (hours)')

%% Reading in and plotting T&H files
for i = sampleNums
cd(sample(i).dir);
sample(i).TH = csvread(strcat('Sample',num2str(i),'Edit.txt'),0,2);%Reads T& H in farenheight and %H from the edited version(that only has the points we need)
sample(i).tempC = ((sample(i).TH(:,1) - 32)/1.8);
sample(i).thAxis = linspace(0,size(sample(i).TH,1)/60,size(sample(i).TH,1));

maxH = max(sample(i).TH(:,2));
maxT = max(sample(i).tempC);
sample(i).HNorm = sample(i).TH(:,2)/maxH(1);
sample(i).TNorm = sample(i).tempC(:)/maxT(1);
% for j = 2:size(dataReal,2);
%     dataNorm = [dataNorm, dataReal(:,j)/maxS(j)];
% end

end




figure(9)
hold on
axis([0 30 42 60])
plot(sample(1).thAxis,sample(1).TH(:,2))% these will just plot the Relative humidities
plot(sample(2).thAxis,sample(2).TH(:,2),'-m')
plot(sample(3).thAxis,sample(3).TH(:,2),'-c')
plot(sample(4).thAxis,sample(4).TH(:,2),'-r')
plot(sample(5).thAxis,sample(5).TH(:,2),'-k')
plot(sample(6).thAxis,sample(6).TH(:,2),'-g')
plot(sample(13).thAxis,sample(13).TH(:,2),'c','LineWidth',1.5,'LineStyle','--')
xlabel('Time(hours)')
ylabel('Percent Humidity')
legend('Sample1 No Additives','Sample2 Ash Only','Sample3 Brick Dust Only','Sample4 Clay Only','Sample5 All Additives','Sample6 Clay/Brick Dust','Sample 13 1% Brick dust')
title('Percent Humidities One to One')

figure(10)
hold on
axis([0 30 42 60])
plot(sample(72).thAxis,sample(72).TH(:,2))% these will just plot the Relative humidities
plot(sample(8).thAxis,sample(8).TH(:,2),'-m')
plot(sample(92).thAxis,sample(92).TH(:,2),'-c')
plot(sample(10).thAxis,sample(10).TH(:,2),'-r')
plot(sample(112).thAxis,sample(112).TH(:,2),'-k')
plot(sample(12).thAxis,sample(12).TH(:,2),'-g')
xlabel('Time(hours)')
ylabel('Percent Humidity')
legend('Sample7 No Additives','Sample8 Ash Only','Sample9 Brick Dust Only','Sample10 Clay Only','Sample11 All Additives','Sample12 Clay/Brick Dust')
title('Percent Humidities Two to One')

figure(11)
hold on
axis([0 30 20 28])
plot(sample(1).thAxis,sample(1).tempC(:))% these will just plot the Temperatures
plot(sample(2).thAxis,sample(2).tempC(:),'-m')
plot(sample(3).thAxis,sample(3).tempC(:),'-c')
plot(sample(4).thAxis,sample(4).tempC(:),'-r')
plot(sample(5).thAxis,sample(5).tempC(:),'-k')
plot(sample(6).thAxis,sample(6).tempC(:),'-g')
plot(sample(13).thAxis,sample(13).tempC(:),'c','LineWidth',1.5,'LineStyle','--')
xlabel('Time(hours)')
ylabel('Degrees Celsius')
legend('Sample1 No Additives','Sample2 Ash Only','Sample3 Brick Dust Only','Sample4 Clay Only','Sample5 All Additives','Sample6 Clay/Brick Dust','Sample 13 1% Brick Dust')
title('Temperature in Degrees C One to One')

figure(12)
hold on
axis([0 30 20 28])
plot(sample(72).thAxis,sample(72).tempC(:))% these will just plot the Temperatures
plot(sample(8).thAxis,sample(8).tempC(:),'-m')
plot(sample(92).thAxis,sample(92).tempC(:),'-c')
plot(sample(10).thAxis,sample(10).tempC(:),'-r')
plot(sample(112).thAxis,sample(112).tempC(:),'-k')
plot(sample(12).thAxis,sample(12).tempC(:),'-g')
xlabel('Time(hours)')
ylabel('Degrees Celsius')
legend('Sample7 No Additives','Sample8 Ash Only','Sample9 Brick Dust Only','Sample10 Clay Only','Sample11 All Additives','Sample12 Clay/Brick Dust')
title('Temperature in Degrees C Two to One')

%% Other stuff

% sampleFitsR = fit(1).resid;
% for i = sampleNums
%     for j = 1:size(sample(i).exp,2)
%         
%     sample = [sampleFitsR, fit(i).beta(4)];
%     end
% end
% % 
% figure(1) % this will plot the normalized amplitudes so that the drying curve can easily be seen
% % subplot(1,2,1);
% hold on
% axis([0 30 0 1.1])
% plot(sample(1).timePoints,sample(1).ampNorm)% these will just plot the exponents
% plot(sample(2).timePoints,sample(2).ampNorm,'-m')
% plot(sample(3).timePoints,sample(3).ampNorm,'-c')
% plot(sample(4).timePoints,sample(4).ampNorm,'-r')
% plot(sample(5).timePoints,sample(5).ampNorm,'-k')
% plot(sample(6).timePoints,sample(6).ampNorm,'-g')
% 
% plot(sample(1).thAxis,sample(1).HNorm)% these will just plot the Relative humidities
% plot(sample(2).thAxis,sample(2).HNorm,'-m')
% plot(sample(3).thAxis,sample(3).HNorm,'-c')
% plot(sample(4).thAxis,sample(4).HNorm,'-r')
% plot(sample(5).thAxis,sample(5).HNorm,'-k')
% plot(sample(6).thAxis,sample(6).HNorm,'-g')
% 
% plot(sample(1).thAxis,sample(1).TNorm)% these will just plot the Relative humidities
% plot(sample(2).thAxis,sample(2).TNorm,'-m')
% plot(sample(3).thAxis,sample(3).TNorm,'-c')
% plot(sample(4).thAxis,sample(4).TNorm,'-r')
% plot(sample(5).thAxis,sample(5).TNorm,'-k')
% plot(sample(6).thAxis,sample(6).TNorm,'-g')
% 
% 
% title('One to One Samples')
% xlabel('Experiment Time (hours)')
% ylabel('Normalized Amplitude')
% legend('Sample1 No Additives','Sample2 Ash Only','Sample3 Brick Dust Only','Sample4 Clay Only','Sample5 All Additives','Sample6 Clay/Brick Dust')%, '2nd Biexponent')
% 
% figure(2)
% % subplot(1,2,2);
% hold on
% axis([0 30 0 1.1])
% plot(sample(72).timePoints,sample(72).ampNorm)% these will just plot the Amplitudes
% plot(sample(8).timePoints,sample(8).ampNorm,'-m')
% plot(sample(9).timePoints,sample(9).ampNorm,'-c')
% plot(sample(10).timePoints,sample(10).ampNorm,'-r')
% plot(sample(11).timePoints,sample(11).ampNorm,'-k')
% plot(sample(12).timePoints,sample(12).ampNorm,'-g')
% 
% plot(sample(72).thAxis,sample(72).HNorm)% these will just plot the Relative humidities
% plot(sample(8).thAxis,sample(8).HNorm,'-m')
% plot(sample(9).thAxis,sample(9).HNorm,'-c')
% plot(sample(10).thAxis,sample(10).HNorm,'-r')
% plot(sample(11).thAxis,sample(11).HNorm,'-k')
% plot(sample(12).thAxis,sample(12).HNorm,'-g')
% 
% plot(sample(72).thAxis,sample(72).TNorm)% these will just plot the Relative humidities
% plot(sample(8).thAxis,sample(8).TNorm,'-m')
% plot(sample(9).thAxis,sample(9).TNorm,'-c')
% plot(sample(10).thAxis,sample(10).TNorm,'-r')
% plot(sample(11).thAxis,sample(11).TNorm,'-k')
% plot(sample(12).thAxis,sample(12).TNorm,'-g')
% 
% title('Two to One Samples')
% xlabel('Experiment Time (hours)')
% ylabel('Normalized Amplitude')
% legend('Sample7 No Additives','Sample8 Ash Only','Sample9 Brick Dust Only','Sample10 Clay Only','Sample11 All Additives','Sample12 Clay/Brick Dust')%, '2nd Biexponent')

%% Otherstuff (cont'd) comparing the reruns to the original values

% i=5;% changing these allows you to easily change the plotted values
% j=52; % your still going to have to change the names though :*(
% 
% close all
% figure(1)%sample 5-52
% subplot(3,1,1)
% plot(sample(i).timePoints,sample(i).ampNorm,'-k')
% hold on
% axis([0 15 0 1.1])
% plot(sample(j).timePoints,sample(j).ampNorm,'-r')
% xlabel('Experiment time [hours]')
% ylabel('Normalized Amp')
% legend('Sample 5 11days old','Sample52 29 days old')
% 
% subplot(3,1,2)
% plot(sample(i).timePoints,sample(i).exp,'-k')
% hold on
% axis([0 15 0 1.1])
% plot(sample(j).timePoints,sample(j).exp,'-r')
% plot(sample(i).timePoints,sample(i).exp2,'-k')
% plot(sample(j).timePoints,sample(j).exp2,'-r')
% xlabel('Experiment time [hours]')
% ylabel('Biexp Amps')
% legend('Sample 5 11days old','Sample52 29 days old')
% 
% subplot(3,1,3)
% errorbar(sample(i).timePoints,sample(i).t2,sample(i).eT2diffl,sample(i).eT2diffu,'-k')
% hold on
% axis([0 15 0 30])
% errorbar(sample(j).timePoints,sample(j).t2,sample(j).eT2diffl,sample(j).eT2diffu,'-r')
% xlabel('Experiment time [hours]')
% ylabel('T2 time [ms]')
% legend('Sample 5 11days old','Sample52 29 days old')
% 
% i=7;% changing these allows you to easily change the plotted values
% j=72; % your still going to have to change the names though :*(
% 
% figure(2)%sample 7-72
% subplot(3,1,1)
% plot(sample(i).timePoints,sample(i).ampNorm,'-k')
% hold on
% axis([0 15 0 1.1])
% plot(sample(j).timePoints,sample(j).ampNorm,'-r')
% xlabel('Experiment time [hours]')
% ylabel('Normalized Amp')
% legend('Sample 7 60.1% H20','Sample72 50.5% H20')
% 
% subplot(3,1,2)
% plot(sample(i).timePoints,sample(i).exp,'-k')
% hold on
% axis([0 15 0 1.1])
% plot(sample(j).timePoints,sample(j).exp,'-r')
% plot(sample(i).timePoints,sample(i).exp2,'-k')
% plot(sample(j).timePoints,sample(j).exp2,'-r')
% xlabel('Experiment time [hours]')
% ylabel('Biexp Amps')
% legend('Sample 7 60.1% H20','Sample72 50.5% H20')
% 
% subplot(3,1,3)
% errorbar(sample(i).timePoints,sample(i).t2,sample(i).eT2diffl,sample(i).eT2diffu,'-k')
% hold on
% axis([0 15 0 30])
% errorbar(sample(j).timePoints,sample(j).t2,sample(j).eT2diffl,sample(j).eT2diffu,'-r')
% xlabel('Experiment time [hours]')
% ylabel('T2 time [ms]')
% legend('Sample 7 60.1% H20','Sample72 50.5% H20')
% 
% i=9;% changing these allows you to easily change the plotted values
% j=92; % your still going to have to change the names though :*(
% 
% figure(3)%sample 9-92
% subplot(3,1,1)
% plot(sample(i).timePoints,sample(i).ampNorm,'-k')
% hold on
% axis([0 15 0 1.1])
% plot(sample(j).timePoints,sample(j).ampNorm,'-r')
% xlabel('Experiment time [hours]')
% ylabel('Normalized Amp')
% legend('Sample 9 46.8% H20','Sample92 49.9% H20')
% 
% subplot(3,1,2)
% plot(sample(i).timePoints,sample(i).exp,'-k')
% hold on
% axis([0 15 0 1.1])
% plot(sample(j).timePoints,sample(j).exp,'-r')
% plot(sample(i).timePoints,sample(i).exp2,'-k')
% plot(sample(j).timePoints,sample(j).exp2,'-r')
% xlabel('Experiment time [hours]')
% ylabel('Biexp Amps')
% legend('Sample 9 46.8% H20','Sample92 49.9% H20')
% 
% subplot(3,1,3)
% errorbar(sample(i).timePoints,sample(i).t2,sample(i).eT2diffl,sample(i).eT2diffu,'-k')
% hold on
% axis([0 15 0 30])
% errorbar(sample(j).timePoints,sample(j).t2,sample(j).eT2diffl,sample(j).eT2diffu,'-r')
% xlabel('Experiment time [hours]')
% ylabel('T2 time [ms]')
% legend('Sample 9 46.8% H20','Sample92 49.9% H20')
% 
% i=11;% changing these allows you to easily change the plotted values
% j=112; % your still going to have to change the names though :*(
% 
% figure(4)%sample 9-92
% subplot(3,1,1)
% plot(sample(i).timePoints,sample(i).ampNorm,'-k')
% hold on
% axis([0 15 0 1.1])
% plot(sample(j).timePoints,sample(j).ampNorm,'-r')
% xlabel('Experiment time [hours]')
% ylabel('Normalized Amp')
% legend('Sample 11 68.2% H20','Sample112 49.7% H20')
% 
% subplot(3,1,2)
% plot(sample(i).timePoints,sample(i).exp,'-k')
% hold on
% axis([0 15 0 1.1])
% plot(sample(j).timePoints,sample(j).exp,'-r')
% plot(sample(i).timePoints,sample(i).exp2,'-k')
% plot(sample(j).timePoints,sample(j).exp2,'-r')
% xlabel('Experiment time [hours]')
% ylabel('Biexp Amps')
% legend('Sample 11 68.2% H20','Sample112 49.7% H20')
% 
% subplot(3,1,3)
% errorbar(sample(i).timePoints,sample(i).t2,sample(i).eT2diffl,sample(i).eT2diffu,'-k')
% hold on
% axis([0 15 0 30])
% errorbar(sample(j).timePoints,sample(j).t2,sample(j).eT2diffl,sample(j).eT2diffu,'-r')
% xlabel('Experiment time [hours]')
% ylabel('T2 time [ms]')
% legend('Sample 11 68.2% H20','Sample112 49.7% H20')