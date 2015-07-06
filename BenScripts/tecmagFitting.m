close all
clear all
addpath(genpath('Z:\TKM\'));

dir1 = ('C:\Users\bmfortman\Documents\Data\MortarDryingTecMag');
dirP25 = ('C:\Users\bmfortman\Documents\Data\Tecmagtest\PinkPetKea25');
dirP5 = ('C:\Users\bmfortman\Documents\Data\Tecmagtest\PinkPet5');
dirH2o = ('C:\Users\bmfortman\Documents\Data\MortarDiffusion');
cd(dirH2o)

%% reading Data from a tecmag file
[ap,s1,s2] = readTecmag('DIWater_n32_nE512_12June2014.tnt'); %reads the data from the file all are complex points
%s1 and s2 are the complex points, s1 is the first entry on s2
sReal = real(s1);
sImag = imag(s1);



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
s2Ac = s2b;% cutting down experiment to actual run amount of 2d experiments
s3 = reshape(s2Ac,64,512); % # acq pts, # echoes, # 2dexperiments
% surf(abs(s3(:,:,3)));
s2sum = sum(s3,1);
s2sum = reshape(s2sum,512,1);%summed data, #Echoes, #2d experiments gives you summed echoes for all echoes in the train
s2sumAbs = abs(s2sum);
% surf(s2sum)
% surf(abs(s2sum))
echoAxis = linspace(150,512*150,512); %gives your time axis #echotime, #number echoes*Echotime, numberEchoes(assuming dwell time of 1u)
% scatter(echoAxis,s2sum(:,1))

% plot(echoAxis(2:end),(s2sum(2:end,1)),'-r')
% hold on% used to check for noise

%% Normalizing data
maxS = max(s2sumAbs);
l = size(s2sum);
sNorm = s2sum(:,1)/maxS(1);
for j = 2:l(2)
    sNorm = [sNorm, s2sum(:,j)/maxS(j)];
end

%% Fitting with Bi or Mono fit
% for sample 7-2 [1500000, 15000]
guesses = [1, 20000]%, 1, 30000];%, 0.1, 7];% remember to change guesses for biexp vs monoexp 
guesses2 = [0,0];%,0,0];
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');
expTime = (ap.ns * 1 * l(2))/3600; %#scans*reptime*#2dScans(taken from size of summed data) /3600 gives hours; /60 gives minutes
timeScale = linspace(0,expTime,l(2)); % 
for j = 1 : l(2) %this is the exp fit t2bifit needs to be changed to t2monofit to convert between the two
fit(j).beta=[0,0]%,0,0];
    while fit(j).beta~=guesses %generates while loop if guess is off
    [fit(j).beta,fit(j).resid,fit(j).J] = nlinfit(echoAxis(3:end)',((sNorm(3:end,j))),@t2monofit_simple,guesses,fitopts);%taking from column 20 for testing
    %catch [fit(j).beta,fit(j).resid,fit(j).J] = nlinfit(echoAxis,((s2sumAbs(:,j))),@t2monofit_simple,guesses2,fitopts);%taking from column 20 for testing
    %end
    guesses=fit(j).beta; %reallocates guess to new result
end
        
    fit(j).pred = t2monofit_simple(fit(j).beta,echoAxis(3:end)');
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

for j = 1:l(2)
% figure(j)
% hold on
% plot(echoAxis',s2sumAbs(:,j))
% plot(echoAxis',fit(j).pred)
% xlabel('useconds')
% ylabel('amplitude')
amplitude(j) = (fit(j).beta(1));

t2Time(j) = (fit(j).beta(2));
end
%% plotting of amplitude/T2's
figure(1)
hold on
plot(timeScale,(amplitude./ap.ns))% dividing by #of scans to get the points for each individual point
xlabel('hours')
ylabel('Amplitude')
legend('Amplitude')

figure(2);
hold on;
% axis([0 14 0 70]) % for sample 7-2
plot(timeScale,(t2Time./ap.ns))% dividing by # of scans to get the points for each individual point
xlabel('Hours')
ylabel('T2time')
legend('T2Time')



