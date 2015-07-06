close all
clear all
addpath(genpath('Z:\TKM\'));

dirP25 = ('C:\Users\bmfortman\Documents\Data\Tecmagtest\PinkPetKea25');
dirP5 = ('C:\Users\bmfortman\Documents\Data\Tecmagtest\PinkPet5');

%% reading Data from a tecmag file
% cd(dirP25); [ap,s1,s2] = readTecmag('tecmagtestPinkPet25.tnt'); % reads the data for the 25
cd(dirP5); [ap,s1,s2] = readTecmag('Pm5TecmagPinkPet.tnt'); % reads the data for the 5
%s1 and s2 are the complex points, s1 is the first entry on s2
sReal = real(s2);
sImag = imag(s2);



% for i=1:80
%     for j=1:100
% intEcho(i,j)=sum(abs(s2(i,j*64-63:j*64)));
%     end
% end

% figure(1)
% hold on
% plot(sReal(1,:)) % plots the first echo set of the data
% scatter(time,intEcho(1,:),'o')
%% summing echoes and reshaping

s2b=sReal';%
s2Ac = s2b(:,(1:12));% cutting down experiment to actual run amount of 2d experiments
s3 = reshape(s2Ac,64,64,12); % # acq pts, # echoes, # 2dexperiments
% surf(abs(s3(:,:,3)));
s2sum = sum(s3,1);
s2sum = reshape(s2sum,64,12);%summed data, #Echoes, #2d experiments gives you summed echoes for all echoes in the train
s2sumAbs = abs(s2sum);
% surf(s2sum)
% surf(abs(s2sum))
echoAxis = linspace(150,64*150,64); %gives your time axis #echotime, #number echoes*Echotime, numberEchoes(assuming dwell time of 1u)
% scatter(echoAxis,s2sum(:,1))

% plot(echoAxis(2:end),(s2sum(2:end,1)),'-r')
% hold on% used to check for noise

%% Normalizing dataTM
maxS = max(s2sumAbs);
l = size(s2sum);
sNorm = s2sum(:,1)/maxS(1);
for j = 2:l(2)
    sNorm = [sNorm, s2sum(:,j)/maxS(j)];
end

%% Fitting with Bi or Mono fit
guesses = [8500, 20000];%, 0.1, 7];% remember to change guesses for biexp vs monoexp 
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');
expTime = (ap.ns * 1 * l(2))/3600; %#scans*reptime*#2dScans(taken from size of summed data) /3600 gives hours; /60 gives minutes
timeScale = linspace(0,expTime,l(2)); % 
for j = 1 : l(2) %this is the exp fit t2bifit needs to be changed to t2monofit to convert between the two
fit(j).beta=[0,0];
    while fit(j).beta~=guesses %generates while loop if guess is off
    [fit(j).beta,fit(j).resid,fit(j).J] = nlinfit(echoAxis',((sNorm(:,j))),@t2monofit_simple,guesses,fitopts);%taking from column 20 for testing
    %catch [fit(j).beta,fit(j).resid,fit(j).J] = nlinfit(echoAxis,((s2sumAbs(:,j))),@t2monofit_simple,guesses2,fitopts);%taking from column 20 for testing
    %end
    guesses=fit(j).beta; %reallocates guess to new result
end
        
    fit(j).pred = t2monofit_simple(fit(j).beta,echoAxis');
%     figure(j)
%     hold on % commented out, plotting of data as it is prepared
%     plot((data1(j).sig(:,1)),(data1(j).sig(:,i+1)))%col 20 for testing
%     plot((data1(j).sig(:,1)),fit(j).pred,'-k')
    %[Ypred,delta] = nlpredci(@t2bifit_simple,echoVector',fit(j).beta,fit(j).resid,'Jacobian',fit(j).J);
    ci = nlparci(fit(j).beta,fit(j).resid,'jacobian',fit(j).J);
%     error(j).Ypreds = (sum(Ypred)/(length(Ypred)));
%     error(j).deltas = (sum(delta)/(length(delta)));
    error(j).ci = ci;
    %error(j).delta = delta;
    %bigFit(j).fit = fit.beta;
    %bigError(j).error = error;

end 



for j = 1:12%l(2)
% figure(j)
% hold on
% plot(echoAxis',s2sumAbs(:,j))
% plot(echoAxis',fit(j).pred)
% xlabel('useconds')
% ylabel('amplitude')
amplitude(j) = (fit(j).beta(1));


t2Time(j) = (fit(j).beta(2));
end

tecmagFits = fit(1).pred; % creates an accumulated fit for the tecmag data
for j = 2:12
tecmagFits = [tecmagFits, fit(j).pred];
end

% figure(1)
% hold on
% plot(timeScale,(amplitude./ap.ns))% dividing by #of scans to get the points for each individual point
% xlabel('hours')
% ylabel('Amplitude')
% legend('Amplitude')
% 
% figure(2);
% hold on;
% % axis([0 14 0 70]) % for sample 7-2
% plot(timeScale,(t2Time./ap.ns))% dividing by # of scans to get the points for each individual point
% xlabel('Hours')
% ylabel('T2time')
% legend('T2Time')


%% Mortar Drying Matlab
tmTest25 = 'C:\Users\bmfortman\Documents\Data\Tecmagtest\PinkPetKea25\1';
tmTest5 = 'C:\Users\bmfortman\Documents\Data\Tecmagtest\PinkPet5\1';

cd(tmTest5); %when changing directories you also need to change the name of the file properly in lines 30,31,38,39

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

% for j = 1:6% For the pm25
%         
%     dataR = strcat('PinkPetKea250',num2str(j),'-decaysRe.dat');
%     dataI = strcat('PinkPetKea250',num2str(j),'-decaysIm.dat');
%     data(j).real = load(dataR);
%     data(j).imag = load(dataI);
% end

for j = 1:6% for the pm5
        
    dataR = strcat('PinkPet50',num2str(j),'-decaysRe.dat');
    dataI = strcat('PinkPet50',num2str(j),'-decaysIm.dat');
    data(j).real = load(dataR);
    data(j).imag = load(dataI);
end

echoAxis = [data(1).real(:,1)]; % works for the kea, not so for the tecMag
dataReal = [data(1).real(:,2), data(1).real(:,3)]; % these are starting at the 4th point to avoid extra noise from the pulse sequence
dataImag = [data(1).imag(:,2), data(1).imag(:,3)]; % the 2nd and 3rd columns are depths, for sample 7 it is all at the same depth so these are disregarded
    
% dataReal= [data(1).real]; % puts all of the data into a single matrix? NOOOOOOOOOOOOOOOOOOOOOOOOOOO
% dataImag= [data(1).imag]; % ALL it does is just start the matrix

for j = 2:length(expNums) % THIS puts all of the data into a single matrix the first column of which is the echoAxis
    dataReal = [dataReal, data(j).real(:,2), data(j).real(:,3)];
    dataImag = [dataImag, data(j).imag(:,2), data(j).imag(:,3)];
end

%% Normalizing data Kea
maxD = max(dataReal);
l = size(s2sum);
dNorm = dataReal(:,1)/maxD(1);
for j = 2:l(2)
    dNorm = [dNorm, dataReal(:,j)/maxD(j)];
end


%% exponential fit

guesses = [1, 2];%,% 0.1, 7]; %these work well for sample 7
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');


for j = 1:size(dNorm,2);% this is to move along the matrices taking the next point
    
    [fit(j).beta,fit(j).resid,fit(j).J] = nlinfit(echoAxis,dNorm(:,j),@t2monofit_simple,guesses,fitopts);% mono exp fit
    fit(j).pred = t2monofit_simple(fit(j).beta,echoAxis); % produces predicted values for easy visualization of fit
%     figure(j) % used for visualizing the fit
%     hold on
%     plot(echoAxis,dataReal(:,j))
%     plot(echoAxis,fit(j).pred,'-k')
    ci = nlparci(fit(j).beta,fit(j).resid,'jacobian',fit(j).J);
    error(j).ci = ci; % gives confidence intervals
  
end

for j = 1:size(dNorm,2); % converts the structures to matrices for ease of graphing
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

keaFits = fit(1).pred; % creates an accumulated fit for the tecmag data
for j = 2:12
keaFits = [keaFits, fit(j).pred];
end

%% plotting Monoexpfit vs. time and depth
expTime = params.repTime * params.nrScans * params.nrExp * params.nrSteps/3600000;% for seconds/1000; for minutes/60000; for hours /3600000
timePoints = linspace(0,expTime,size(exp,2));

% figure(1)
% hold on
% plot(timePoints,exp)%,eExpdiffl,eExpdiffu)
% %plot(expDepth,biExp2,'-r')
% xlabel('Experiment Time (hours)')
% ylabel('Amplitude')
% legend('1st Exponent')%, '2nd Biexponent')
% 
% figure(2)
% hold on
% plot(timePoints,t2)%,eT2diffl,eT2diffu)
% %plot(expDepth,biExpt22,'-r')
% xlabel('Experiment Time (hours)')
% ylabel('T2 time')
% legend('1st T2 time')%, '2nd T2 time')

%% THE DIRECT COMPARISON

for j = 1:12
    figure(j)
    hold on
    plot(echoAxis,sNorm(:,j))
    plot(echoAxis,dNorm(:,j),'-k')
    plot(echoAxis,tecmagFits(:,j))
    plot(echoAxis,keaFits(:,j),'-k')
    legend('Tecmag','Kea')
end
dataBoth = [sNorm(:,1),dNorm(:,1)];
for j = 2:12
    figure(13)
    dataBoth = [dataBoth, sNorm(:,j), dNorm(:,j)];
end

surf(dataBoth)
title('Tecmag starts on one, then alternates with kea')
       
