close all
clear
clc
%% setting up parameters

nEchoes = 512
tE = 150e-6
nrPts = 69
nrBlank = 5 
%% setting up directories

% dir1= ('C:\Users\bmfortman\Documents\Data\MortarCuring\Brick Dust Experiments');% for the lab computer
dir1 = ('C:\Users\benjamin\Documents\Data\Mortar Curing\Brick Dust Experiments'); % personal computer

dir2 = ('C:\Users\benjamin\Documents\Data\Mortar Curing\Post cure Brick dust Exp');% personal computer directory for cured mortar CPMG's

% there is another dir 2 below that will need to be changed
cd(dir1)
%% reading in data

% setting up the names for the files
sample(1).name = ('Sample0_10292014_256scans_512Echoes_CPMG.tnt');
sample(2).name = ('Sample1_10212014_256scans_512Echoes_CPMG.tnt');
sample(3).name = ('Sample2_10222014_256scans_512Echoes_CPMG.tnt');
sample(4).name = ('Sample3_10232014_256scans_512Echoes_CPMG.tnt');
sample(5).name = ('Sample4_10242014_256scans_512Echoes_CPMG.tnt');
sample(6).name = ('Sample5_10252014_256scans_512Echoes_CPMG.tnt');
sample(7).name = ('Sample6_10262014_256scans_512Echoes_CPMG.tnt');
sample(8).name = ('Sample7_10272014_256scans_512Echoes_CPMG.tnt');
sample(9).name = ('Sample8_10282014_256scans_512Echoes_CPMG.tnt');


%reading data in, reshaping the data and cutting out blank points
for i = 1:length(sample)
    [sample(i).params,sample(i).data] = readTecmag4d(sample(i).name);
    sample(i).data = reshape(real(sample(i).data),nrPts,nEchoes);
    sample(i).data = sample(i).data(1:end-nrBlank,:);
    sample(i).sumData = sum(sample(i).data);
    sample(i).normSum = sample(i).sumData(:)./max(sample(i).sumData);
end


%% fitting with nlinfit and such
echoAxis = (1:nEchoes)*tE;
guesses = [1,0.03,0.1,0.001];

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

%% reading in, then calculation of SV ratio assuming spherical pore
for i = 1:length(sample) % xscalars to convert T2's to pore sizes
    Xscales(i) = load(strcat('BrickdustSample',num2str(i-1),'xScalar.txt'));
end

poreSize = Xscales'.*T2; % gives values for the different pore sizes
smallPoreSize = Xscales'.*T22; % gives values for the smaller pore size
rhoVals = Xscales./(1e6*3*2); % converts X scalar values back into the rho vals(rho = 1/Harmonic mean)* sVratio^-1)

%% fitting to the single fixed, biexponential fxn
guesses2 = [0.80,0.035];
for i = 1:length(sample)
    [fixed(i).beta,fixed(i).resid,fixed(i).j] = nlinfit(echoAxis(3:end),sample(i).normSum(3:end)',@t2bifit_monofixed,guesses2);
    fixed(i).pred = t2bifit_monofixed(fixed(i).beta,echoAxis(3:end));
    fixed(i).ci = nlparci(fixed(i).beta,fixed(i).resid,'jacobian',fixed(i).j);

end

%again compiling fits
FT2 = fixed(1).beta(2);
Famp = fixed(1).beta(1);
FT2eu = fixed(i).ci(2,:);
for i = 2:length(sample)
    FT2 = [FT2;fixed(i).beta(2)];
    Famp = [Famp; fixed(i).beta(1)];
    FT2eu = [FT2eu; fixed(i).ci(2,:)];
end

%error stuff to get it in the right format
FT2eu = abs([FT2eu(:,1) - FT2,FT2eu(:,2) - FT2]);

%% Plotting all of the functions

dustPer = ((1:length(sample))-1)*0.25; %creates per of brick dust

figure(1)
plot(dustPer,amp,'r')
hold on
plot(dustPer,amp2,'k')
plot(dustPer,Famp,'b')
xlabel('percentage Brick dust')
ylabel('Normalized amp')
legend('1st amplitude','2nd amplitude','amplitude fitted to mean of 2nd Amplitude')

figure(2)
errorbar(dustPer,T2,T2eu(:,1),T2eu(:,2),'r')
hold on
errorbar(dustPer,T22,T22eu(:,1),T22eu(:,2),'k')
errorbar(dustPer,FT2,FT2eu(:,1),FT2eu(:,2),'b')
xlabel('percentage Brick dust')
ylabel('Normalized T_2 [s]')
legend('1st T_2','2nd T_2','T_2 fitted to mean of 2nd T_2')

figure(3)
plot(dustPer,poreSize,'b')
axis([0 2 0 22])
hold on
plot(dustPer,smallPoreSize,'r')
ylabel('pore diameter [um]')
xlabel('BrickDust Percentage')
legend('Larger Pore size','Smaller pore size')
title('Average Pore size based upon T2 Values')
%% Clearing the data before moving on

% clear
% dir2 = ('C:\Users\benjamin\Documents\Data\Mortar Curing\Post cure Brick dust Exp');% personal computer directory for cured mortar CPMG's


%% Move onto the second set of post cure data
cd(dir2)

sample(1).name = ('Sample0_1222014_256scans_512Echoes_CPMG.tnt');
sample(2).name = ('Sample1_12092014_256scans_512Echoes_CPMG.tnt');
sample(3).name = ('Sample2CPMG12102014.tnt');
sample(4).name = ('Sample3CPMG mortar12112014.tnt');
sample(5).name = ('Sample4CPMG mortar12122014.tnt');
sample(6).name = ('Sample5_CPMG_12132014.tnt');
sample(7).name = ('Sample6_CPMG_12142014.tnt');
sample(8).name = ('Sample7_CPMG_12152014.tnt');
sample(9).name = ('Sample8_CPMG_12162014.tnt');
%% other sample names
% cd('C:\Users\benjamin\Documents\Data\Mortar Curing\Brick Dust Experiments\Retest CPMGs')
% sample(1).name = ('Sample02_CPMG_242015.tnt');
% sample(2).name = ('Sample12_CPMG_242015.tnt');
% sample(3).name = ('Sample22_CPMG_242015.tnt');
% sample(4).name = ('Sample32_CPMG_242015.tnt');
% sample(5).name = ('Sample42_CPMG_242015.tnt');
% sample(6).name = ('Sample52_CPMG_242015.tnt');
% sample(7).name = ('Sample62_CPMG_242015.tnt');
% sample(8).name = ('Sample72_CPMG_242015.tnt');
% sample(9).name = ('Sample82_CPMG_242015.tnt');

% sample(1).name = ('Sample03_CPMG_252015.tnt');
% sample(2).name = ('Sample13_CPMG_252015.tnt');
% sample(3).name = ('Sample23_CPMG_252015.tnt');
% sample(4).name = ('Sample33_CPMG_252015.tnt');
% sample(5).name = ('Sample43_CPMG_252015.tnt');
% sample(6).name = ('Sample53_CPMG_252015.tnt');
% sample(7).name = ('Sample63_CPMG_252015.tnt');
% sample(8).name = ('Sample73_CPMG_242015.tnt');
% sample(9).name = ('Sample83_CPMG_242015.tnt');

%reretests
% clear('sample')
% sample(1).name = ('Sample5_CPMG_2112015-1.tnt');
% sample(2).name = ('Sample1_CPMG_2112015.tnt');
% sample(3).name = ('Sample1_CPMG_2112015-2.tnt');
% sample(4).name = ('Sample1_CPMG_2112015-3.tnt');
% sample(5).name = ('Sample1_CPMG_2132015-1.tnt');
% sample(6).name = ('Sample1_CPMG_2132015-2.tnt');
% sample(7).name = ('Sample1_CPMG_2132015-3.tnt');
% sample(8).name = ('Sample5_CPMG_2132015-1.tnt');
% sample(9).name = ('Sample5_CPMG_2132015-2.tnt');
% sample(10).name = ('Sample5_CPMG_2132015-3.tnt');


%reading data in, reshaping the data and cutting out blank points
for i = 1:length(sample)
    [sample(i).params,sample(i).data] = readTecmag4d(sample(i).name);
    sample(i).data = reshape(real(sample(i).data),nrPts,nEchoes);
    sample(i).data = sample(i).data(1:end-nrBlank,:);
    sample(i).sumData = sum(sample(i).data);
    sample(i).normSum = sample(i).sumData(:)./max(sample(i).sumData);
end


%% fitting with nlinfit and such
echoAxis = (1:nEchoes)*tE;
guesses = [1,0.03,0.1,0.001];

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

%% fitting to the single fixed, biexponential fxn
guesses2 = [0.80,0.035];
for i = 1:length(sample)
    [fixed(i).beta,fixed(i).resid,fixed(i).j] = nlinfit(echoAxis(3:end),sample(i).normSum(3:end)',@t2bifit_monofixed,guesses2);
    fixed(i).pred = t2bifit_monofixed(fixed(i).beta,echoAxis(3:end));
    fixed(i).ci = nlparci(fixed(i).beta,fixed(i).resid,'jacobian',fixed(i).j);

end

%again compiling fits
FT2 = fixed(1).beta(2);
Famp = fixed(1).beta(1);
FT2eu = fixed(i).ci(2,:);
for i = 2:length(sample)
    FT2 = [FT2;fixed(i).beta(2)];
    Famp = [Famp; fixed(i).beta(1)];
    FT2eu = [FT2eu; fixed(i).ci(2,:)];
end

%error stuff to get it in the right format
FT2eu = abs([FT2eu(:,1) - FT2,FT2eu(:,2) - FT2]);

%% Plotting all of the functions

dustPer = ((1:length(sample))-1)*0.25; %creates per of brick dust

figure(1)
plot(dustPer,amp,'r','LineWidth',2,'LineStyle','-.')
hold on
plot(dustPer,amp2,'k','LineWidth',2,'LineStyle','-.')
plot(dustPer,Famp,'b','LineWidth',2,'LineStyle','-.')
xlabel('percentage Brick dust')
ylabel('Normalized amp')
legend('1st amplitude','2nd amplitude','amplitude fitted to mean of 2nd Amplitude','1st amplitude Postcure','2nd amplitude Postcure','amplitude fitted to mean of 2nd Amplitude Postcure')

figure(2)
errorbar(dustPer,T2,T2eu(:,1),T2eu(:,2),'r','LineWidth',2,'LineStyle','-.')
hold on
errorbar(dustPer,T22,T22eu(:,1),T22eu(:,2),'k','LineWidth',2,'LineStyle','-.')
errorbar(dustPer,FT2,FT2eu(:,1),FT2eu(:,2),'b','LineWidth',2,'LineStyle','-.')
xlabel('percentage Brick dust')
ylabel('Normalized T_2 [s]')
legend('1st T_2','2nd T_2','T_2 fitted to mean of 2nd T_2','1st T_2Post Cure','2nd T_2 Postcure','T_2 fitted to mean of 2nd T_2 Postcure')