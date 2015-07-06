clear
close all
clc
addpath(genpath('Z:\TKM\'));

%% Mortar Drying Matlab
maindir = 'C:\Users\bmfortman\Documents\Data\';

dir11 = 'C:\Users\bmfortman\Documents\Data\MortarDrying\OneToOne\1\';
dir21 = 'C:\Users\bmfortman\Documents\Data\MortarDrying\TwoToOne\1\';
sample1 = 'C:\Users\bmfortman\Documents\Data\MortarDrying\1-1Samp1\1'
sample7= 'C:\Users\bmfortman\Documents\Data\MortarDrying\2-1Samp7\1';
tmTest25 = 'C:\Users\bmfortman\Documents\Data\Tecmagtest\PinkPetKea25\1';
tmTest5 = 'C:\Users\bmfortman\Documents\Data\Tecmagtest\PinkPet5\1';

cd(tmTest25); %when changing directories you also need to change the name of the file properly in lines 30,31,38,39

%% Reading in Paramaters

params.echoTime = readpar_Kea('acqu.par','echoTime');
params.pulseLength = readpar_Kea('acqu.par','pulseLength');
params.dwellTime = readpar_Kea('acqu.par','dwellTime');
params.nrExp = readpar_Kea('acqu.par','nrExp');
params.nrSteps = readpar_Kea('acqu.par','nrSteps');
params.nrScans = readpar_Kea('acqu.par','nrScans');
params.repTime = readpar_Kea('acqu.par','repTime');

%% reading mortar data in then plotting
expNums = [1 : (params.nrExp)]; % should cover the number of experiments

for j = 1:6%9% works provided there are more than 9 experiments
        
    dataR = strcat('PinkPetKea250',num2str(j),'-decaysRe.dat');
    dataI = strcat('PinkPetKea250',num2str(j),'-decaysIm.dat');
    data(j).real = load(dataR);
    data(j).imag = load(dataI);
end
% 
% for j = 10:length(expNums) % covers the rest of the experiments
%         
%     dataR = strcat('2-1Samp7',num2str(j),'-decaysRe.dat'); % need to change to name of one of the directory file for different experiments
%     dataI = strcat('2-1Samp7',num2str(j),'-decaysIm.dat');
%     data(j).real = load(dataR);
%     data(j).imag = load(dataI);
% end

echoAxis = [data(1).real(:,1)]; % works for the kea, not so for the tecMag
dataReal = [data(1).real(:,2), data(1).real(:,3)]; % these are starting at the 4th point to avoid extra noise from the pulse sequence
dataImag = [data(1).imag(:,2), data(1).imag(:,3)]; % the 2nd and 3rd columns are depths, for sample 7 it is all at the same depth so these are disregarded
    
% dataReal= [data(1).real]; % puts all of the data into a single matrix? NOOOOOOOOOOOOOOOOOOOOOOOOOOO
% dataImag= [data(1).imag]; % ALL it does is just start the matrix

for j = 2:length(expNums) % THIS puts all of the data into a single matrix the first column of which is the echoAxis
    dataReal = [dataReal, data(j).real(:,2), data(j).real(:,3)];
    dataImag = [dataImag, data(j).imag(:,2), data(j).imag(:,3)];
end

%% exponential fit

guesses = [0.1, 16];%,% 0.1, 7]; %these work well for sample 7
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');


for j = 1:size(dataReal,2);% this is to move along the matrices taking the next point
    
    [fit(j).beta,fit(j).resid,fit(j).J] = nlinfit(echoAxis,dataReal(:,j),@t2monofit_simple,guesses,fitopts);% mono exp fit
    fit(j).pred = t2monofit_simple(fit(j).beta,echoAxis); % produces predicted values for easy visualization of fit
%     figure(j) % used for visualizing the fit
%     hold on
%     plot(echoAxis,dataReal(:,j))
%     plot(echoAxis,fit(j).pred,'-k')
    ci = nlparci(fit(j).beta,fit(j).resid,'jacobian',fit(j).J);
    error(j).ci = ci; % gives confidence intervals
  
end

for j = 1:size(dataReal,2); % converts the structures to matrices for ease of graphing
    exp(j) = fit(j).beta(1);
    t2(j) = fit(j).beta(2);
    eExp(1,j) = error(j).ci(1);
    eExp(2,j) = error(j).ci(3);
    eT2(1,j) = error(j).ci(2);
    eT2(2,j) = error(j).ci(4);
end

eExpdiffl =  abs(eExp(1,:)-exp); %converts confidence intervals into distance from values for use with the errorbar function
eExpdiffu =  (eExp(2,:)-exp);
eT2diffl =  abs(eT2(1,:)-t2);
eT2diffu =  (eT2(1,:)-t2);

%% plotting Monoexpfit vs. time and depth
expTime = params.repTime * params.nrScans * params.nrExp * params.nrSteps/3600000;% for seconds/1000; for minutes/60000; for hours /3600000
timePoints = linspace(0,expTime,size(exp,2));

figure(1)
hold on
plot(timePoints,exp)%,eExpdiffl,eExpdiffu)
%plot(expDepth,biExp2,'-r')
xlabel('Experiment Time (hours)')
ylabel('Amplitude')
legend('1st Exponent')%, '2nd Biexponent')

figure(2)
hold on
plot(timePoints,t2)%,eT2diffl,eT2diffu)
%plot(expDepth,biExpt22,'-r')
xlabel('Experiment Time (hours)')
ylabel('T2 time')
legend('1st T2 time')%, '2nd T2 time')

