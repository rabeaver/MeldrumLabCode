close all
clear
clc
%% setting up parameters

nEchoes = 512
tE = 150e-6
nrPts = 69
nrBlank = 5 
DataperPoint = 3; % the number of files per point
%% setting up directories

% dir1= ('C:\Users\bmfortman\Documents\Data\MortarCuring\Brick Dust Experiments');% for the lab computer
dir1 = ('C:\Users\benjamin\Documents\Data\Mortar Curing\Pozzolon Experiments\CPMG and T2d Postcure'); % personal computer
dir2 = ('C:\Users\benjamin\Documents\Data\Mortar Curing\Pozzolon Experiments\')% for the precure experiments
cd(dir1)
%% reading in data

% setting up the names for the files
data(1).name = ('Sample0Pozz_CPMG_2232015.tnt');
data(2).name = ('Sample0Pozz-2_CPMG_2232015.tnt');
data(3).name = ('Sample0Pozz-3_CPMG_2232015.tnt');
data(4).name = ('Sample1Pozz_CPMG_2242015.tnt');
data(5).name = ('Sample1Pozz-2_CPMG_2242015.tnt');
data(6).name = ('Sample1Pozz-3_CPMG_2242015.tnt');
data(7).name = ('Sample2Pozz_CPMG_2252015.tnt');
data(8).name = ('Sample2Pozz-2_CPMG_2252015.tnt');
data(9).name = ('Sample2Pozz-3_CPMG_2252015.tnt');
data(10).name = ('Sample3Pozz_CPMG_2262015.tnt');
data(11).name = ('Sample3Pozz-2_CPMG_2262015.tnt');
data(12).name = ('Sample3Pozz-3_CPMG_2262015.tnt');
data(13).name = ('Sample4Pozz-1_CPMG_2272015.tnt');
data(14).name = ('Sample4Pozz-2_CPMG_2272015.tnt');
data(15).name = ('Sample4Pozz-3_CPMG_2272015.tnt');
data(16).name = ('Sample5Pozz-1_CPMG_2272015.tnt');
data(17).name = ('Sample5Pozz-2_CPMG_312015.tnt');
data(18).name = ('Sample5Pozz-3_CPMG_312015.tnt');
data(19).name = ('Sample6Pozz_CPMG_322015.tnt');
data(20).name = ('Sample6Pozz-2_CPMG_322015.tnt');
data(21).name = ('Sample6Pozz-3_CPMG_322015.tnt');
data(22).name = ('Sample7Pozz_CPMG_332015.tnt');
data(23).name = ('Sample7Pozz-2_CPMG_332015.tnt');
data(24).name = ('Sample7Pozz-3_CPMG_332015.tnt');
data(25).name = ('Sample8Pozz_CPMG_342015.tnt');
data(26).name = ('Sample8Pozz-2_CPMG_342015.tnt');
data(27).name = ('Sample8Pozz-3_CPMG_342015.tnt');
% the above samples are for the pozzolons post cure, below post cure


%% reading data in, reshaping the data and cutting out blank points
for i = 1:length(data)
    [data(i).params,data(i).data] = readTecmag4d(data(i).name);
    data(i).data = reshape(real(data(i).data),nrPts,nEchoes);
    data(i).data = data(i).data(1:end-nrBlank,:);
end

% taking an average value from every 3 data files, to use later for the
% fitting
for i=1:length(data)/DataperPoint% for every 3datapoints
    sample(i).data= (data(i).data + data(i+1).data + data(i+2).data);
    sample(i).sumData = sum(sample(i).data);
    sample(i).normSum = sample(i).sumData(:)./max(sample(i).sumData);
end

%% fitting with nlinfit and such
echoAxis = (1:nEchoes)*tE;
guesses = [0.8,0.04,0.2,0.005];

for i = 1:length(sample)
    [sample(i).beta,sample(i).resid,sample(i).j] = nlinfit(echoAxis(3:end),sample(i).normSum(3:end)',@t2bifit_simple,guesses);
    sample(i).pred = t2bifit_simple(sample(i).beta,echoAxis(3:end));
    sample(i).ci = nlparci(sample(i).beta,sample(i).resid,'jacobian',sample(i).j);
end

%for compiling fits and finding the average
T2 = sample(1).beta(2);
T22 = sample(1).beta(4);
amp = sample(1).beta(1);
amp2 = sample(1).beta(3);
T2eu = sample(1).ci(2,:);
T22eu = sample(1).ci(4,:);

for i = 2:length(sample)
    T2 = [T2;sample(i).beta(2)];
    T22 = [T22;sample(i).beta(4)];
    amp = [amp; sample(i).beta(1)];
    amp2 = [amp2; sample(i).beta(3)];
    T2eu = [T2eu; sample(i).ci(2,:)];
    T22eu = [T22eu; sample(i).ci(4,:)];

end

%error stuff to get it in the right format
T2eu = abs([T2eu(:,1) - T2,T2eu(:,2) - T2]);
T22eu = abs([T22eu(:,1) - T22,T22eu(:,2) - T22]);

avgT2 = mean(T2);
avgT22 = mean(T22);
avgAmp = mean(amp);
avgAmp2 = mean(amp2);
%% Plotting all of the functions

dustPer = ((1:length(sample))-1)*0.25; %creates per of brick dust

figure(1)
plot(dustPer,amp,'r')
hold on
plot(dustPer,amp2,'k')
xlabel('percentage Brick dust')
ylabel('Normalized amp')
legend('1st amplitude','2nd amplitude')

figure(2)
errorbar(dustPer,T2,T2eu(:,1),T2eu(:,2),'r')
hold on
errorbar(dustPer,T22,T22eu(:,1),T22eu(:,2),'k')
xlabel('percentage Brick dust')
ylabel('Normalized T_2 [s]')
legend('1st T_2','2nd T_2')

%% Other processing of other cpmgs
clear('data')
cd(dir2)
data(1).name = ('Sample1Pozzolons_CPMG_1302015.tnt');
data(2).name = ('Sample4Pozzolons_CPMG_1302015.tnt');
data(3).name = ('Sample8Pozzolons_CPMG_212015.tnt');

for i = 1:length(data)
    [data(i).params,data(i).data] = readTecmag4d(data(i).name);
    data(i).data = reshape(real(data(i).data),nrPts,nEchoes);
    data(i).data = data(i).data(1:end-nrBlank,:);
    data(i).sumData = sum(data(i).data);
    data(i).normSum = data(i).sumData(:)./max(data(i).sumData);
end


for i = 1:length(data)
    [data(i).beta,data(i).resid,data(i).j] = nlinfit(echoAxis(3:end),data(i).normSum(3:end)',@t2bifit_simple,guesses);
    data(i).pred = t2bifit_simple(data(i).beta,echoAxis(3:end));
    data(i).ci = nlparci(data(i).beta,data(i).resid,'jacobian',data(i).j);
end

%for compiling fits and finding the average
T2 = data(1).beta(2);
T22 = data(1).beta(4);
amp = data(1).beta(1);
amp2 = data(1).beta(3);
T2eu = data(1).ci(2,:);
T22eu = data(1).ci(4,:);

for i = 2:length(data)
    T2 = [T2;data(i).beta(2)];
    T22 = [T22;data(i).beta(4)];
    amp = [amp; data(i).beta(1)];
    amp2 = [amp2; data(i).beta(3)];
    T2eu = [T2eu; data(i).ci(2,:)];
    T22eu = [T22eu; data(i).ci(4,:)];

end

%error stuff to get it in the right format
T2eu = abs([T2eu(:,1) - T2,T2eu(:,2) - T2]);
T22eu = abs([T22eu(:,1) - T22,T22eu(:,2) - T22]);

dustPer = [0,1,2];

figure(1)
plot(dustPer,amp,'g')
plot(dustPer,amp2,'g')
legend('1st amplitude','2nd amplitude','Precure amp')

figure(2)
errorbar(dustPer,T2,T2eu(:,1),T2eu(:,2),'g')
errorbar(dustPer,T22,T22eu(:,1),T22eu(:,2),'g')
legend('1st T_2','2nd T_2','precureT_2')
