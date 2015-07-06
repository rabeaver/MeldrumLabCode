clear
close all
clc
addpath(genpath('Z:\TKM\'));

maindir = 'C:\Users\bmfortman\Documents\Data\ArmoryCPMG\';
dirnums = [5 22]
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

%% fitting for a fit

guesses = [2, 4000]%, 1, 1420];
guesses2 = [8, 10163]%, 1, 3000]
timeAxis = linspace(0,8960,256);

fitopts = statset('MaxIter',500,'TolX',1e-14,'UseParallel',true,'Display','off');

figure(1)
hold on
plot(timeAxis,finalData(1).real(:,19))
plot(timeAxis,finalData(2).real(:,17),'-k')


[fit(1).beta,fit(1).resid,fit(1).J] = nlinfit(timeAxis',finalData(1).real(:,19),@t2monofit_simple,guesses,fitopts);
    fit(1).pred = t2monofit_simple(fit(1).beta,timeAxis');
    
    [fit(2).beta,fit(2).resid,fit(2).J] = nlinfit(timeAxis',finalData(2).real(:,17),@t2monofit_simple,guesses2,fitopts);
    fit(2).pred = t2monofit_simple(fit(2).beta,timeAxis');
        ci = nlparci(fit(1).beta,fit(1).resid,'jacobian',fit(1).J);
%     error(j).Ypreds = (sum(Ypred)/(length(Ypred)));
%     error(j).deltas = (sum(delta)/(length(delta)));
    error(1).ci = ci;
    
     ci = nlparci(fit(2).beta,fit(2).resid,'jacobian',fit(2).J);
%     error(j).Ypreds = (sum(Ypred)/(length(Ypred)));
%     error(j).deltas = (sum(delta)/(length(delta)));
    error(2).ci = ci;
    
plot(timeAxis,fit(1).pred)
plot(timeAxis,fit(2).pred,'-k')
plot(timeAxis,fit(2).pred2,'-r')
display (fit(1).beta)
display (fit(2).beta)

