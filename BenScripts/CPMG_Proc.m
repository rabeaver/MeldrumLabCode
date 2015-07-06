clear
close all
clc
addpath(genpath('Z:\TKM\'));

%%
maindir = 'C:\Users\bmfortman\Documents\Data\ArmoryCPMG\';
fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');
%dirnums = [6 7 9 10 11 13 14 15 18 19 20]; %FOR JAMESTOWN i took out day 1 bc of high error
%time = [12 14 16 30 33 40 42 45 47 50 52]; %for Jamestown as i did with its time 1
dirnums = [5 22] %[2 5 6 8 9 10 11 13 17 19 21 22]; %FOR ARMORY I took out 15 and 14 bc of high error values also 20 is having issues with echopoints....
time = [1 10 12 14 16 30 33 35 43 45 50 52]; %for armory and their corresponding times, 40,39
for j = 1:length(dirnums)
    dir = strcat(maindir,num2str(dirnums(j)),'\');
    cd(dir);

data1 = load('data.csv');
data2 = load('data2.csv');

for i = 1:size(data2,2)/2
    data2Re(:,i) = data2(:,i*2-1);
    data2Im(:,i) = data2(:,i*2);
end

dataCp = complex(data2Re,data2Im);
dataCp_phased = autophase(dataCp,1);

params(j).echoTime = readpar_Kea('acqu.par','echoTime');
params(j).pulseLength = readpar_Kea('acqu.par','pulseLength');
params(j).dwellTime = readpar_Kea('acqu.par','dwellTime');
params(j).nrEchoes = readpar_Kea('acqu.par','nrEchoes');
params(j).nrPnts = readpar_Kea('acqu.par','nrPnts');
params(j).nrScans = readpar_Kea('acqu.par','nrScans');
params(j).b1Freq = readpar_Kea('acqu.par','b1Freq');

finalData(j).real = real(dataCp_phased);
finalData(j).imag = imag(dataCp_phased);

    cd(maindir);
    
    clear('data1','data2','data2Re','data2Im','dataCp','dataCp_phased')
end

%% Gaussian filtering
for j = 1:length(dirnums)
x = params(j).dwellTime: params(j).dwellTime: params(j).nrPnts * params(j).dwellTime;
b = 8;
c = params(j).pulseLength/(2*sqrt(2*log(2)));
y = exp(-(x-b).^2./(2*c^2));
plot(x,y)

data(j).longEcho = reshape(finalData(j).real',1,size(finalData(j).real,1)*size(finalData(j).real,2));
Gfilter = repmat(y,params(j).nrEchoes);
Gfilter = Gfilter(1,:);
data(j).filtered = data(j).longEcho.*Gfilter;

clear ('x','Gfilter','y')
plot(data(j).filtered)
data(j).sepEchoes = reshape(data(j).filtered,params(j).nrPnts,params(j).nrEchoes);
data(j).sepEchoes = data(j).sepEchoes';

 figure(j)
 surf(data(j).sepEchoes)
 shading flat

data(j).echoInt = sum(data(j).sepEchoes, 2);
plot(data(j).echoInt);

end

%% biexponential fitting
close all
%guesses = [21, 14000, 12, 3000]; %for jamestown
guesses = [12.9, 14000, 9, 3492];%for armory

for j=1:length(dirnums)
echoVector = 0:params(j).echoTime:(params(j).nrEchoes-1)*params(j).echoTime;

    [fit(j).beta,fit(j).resid,fit(j).J] = nlinfit(echoVector',data(j).echoInt,@t2bifit_simple,guesses,fitopts);
    fit(j).pred = t2bifit_simple(fit(j).beta,echoVector);
    figure(j)
    hold on
    plot(echoVector,data(j).echoInt)
    plot(echoVector,fit(j).pred,'-k')
    %[Ypred,delta] = nlpredci(@t2bifit_simple,echoVector',fit(j).beta,fit(j).resid,'Jacobian',fit(j).J);
    ci = nlparci(fit(j).beta,fit(j).resid,'jacobian',fit(j).J);
%     error(j).Ypreds = (sum(Ypred)/(length(Ypred)));
%     error(j).deltas = (sum(delta)/(length(delta)));
    error(j).ci = ci;
    %error(j).delta = delta;
end


%% Plotting amp v. time/ T2 V time
close all
for j = 1:length(dirnums)
    exp1(j) = fit(j).beta(1);

    exp2(j) = fit(j).beta(3);

    t21(j) = fit(j).beta(2);

    t22(j) = fit(j).beta(4);
    
    errorExp1l(j) = error(j).ci(1);
    errorExp1u(j) = error(j).ci(1,2);
    errorT21l(j) = error(j).ci(2);
    errorT21u(j) = error(j).ci(2,2);
    errorExp2l(j) = error(j).ci(3);
    errorExp2u(j) = error(j).ci(3,2);
    errorT22l(j) = error(j).ci(4);
    errorT22u(j) = error(j).ci(4,2);
    
    
    %errorsd(j) = error(j).deltas;
end

figure(1)
hold on
errorbar(time,exp1,errorExp1l,errorExp1u,'-k')
errorbar(time,exp2,errorExp2l,errorExp2u,'-r')
legend('exp1','exp2')
title('Exponential amplitude comparison')

figure(2)
hold on
errorbar(time,t21,errorT21l,errorT21u,'-k')
errorbar(time,t22,errorT22l,errorT22u,'-r')
legend('1st T2 time', '2nd T2 time')
title('T2 time comparison')

%%other stuff
%nlpredci
