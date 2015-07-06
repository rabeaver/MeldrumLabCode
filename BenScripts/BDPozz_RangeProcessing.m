close all
clear
clc
%% setting up parameters

nEchoes = 512
tE = 150e-6
nrPts = 69
nrBlank = 5 
pozzDust = [0 0.5 1 2 3 4 5 10];% pozzolon in samples
brickDust = [0 1 2 3 4 6 8 10 15 20];% brick dust in samples
% brickDust = [0 1 2 3 4 6 8 10 15 20 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2.0]
%% setting up directories

% dir1= ('C:\Users\bmfortman\Documents\Data\MortarCuring\Brick Dust Experiments');% for the lab computer
dir1 = ('C:\Users\benjamin\Documents\Data\Mortar Curing\B-P Mortar range experiments\'); % personal computer
cd(dir1)
%% reading in data

% Pozzolons are listed first
data(1).name = ('Sample0P(1)-1_CPMG_342015.tnt');
data(2).name = ('Sample0P(1)-2_CPMG_342015.tnt');
data(3).name = ('Sample0.5P(1)-1_CPMG_342015.tnt');
data(4).name = ('Sample0.5P(1)-2_CPMG_342015.tnt');
data(5).name = ('Sample1P(1)-1_CPMG_342015.tnt');
data(6).name = ('Sample1P(1)-2_CPMG_342015.tnt');
data(7).name = ('Sample2P(1)-1_CPMG_342015.tnt');
data(8).name = ('Sample2P(1)-2_CPMG_342015.tnt');
data(9).name = ('Sample3P(1)-1_CPMG_342015.tnt');
data(10).name = ('Sample3P(1)-2_CPMG_342015.tnt');
data(11).name = ('Sample4P(1)-1_CPMG_342015.tnt');
data(12).name = ('Sample4P(1)-2_CPMG_342015.tnt');
data(13).name = ('Sample5P(1)-1_CPMG_342015.tnt');
data(14).name = ('Sample5P(1)-2_CPMG_342015.tnt');
data(15).name = ('Sample10P(1)-1_CPMG_342015.tnt');
data(16).name = ('Sample10P(1)-2_CPMG_342015.tnt');

% Brick Dust is next
data(17).name = ('Sample0B(1)-1_CPMG_352015.tnt');
data(18).name = ('Sample0B(1)-2_CPMG_352015.tnt');
data(19).name = ('Sample1B(1)-1_CPMG_352015.tnt');
data(20).name = ('Sample1B(1)-2_CPMG_352015.tnt');
data(21).name = ('Sample2B(1)-1_CPMG_352015.tnt');
data(22).name = ('Sample2B(1)-2_CPMG_352015.tnt');
data(23).name = ('Sample3B(1)-1_CPMG_352015.tnt');
data(24).name = ('Sample3B(1)-2_CPMG_352015.tnt');
data(25).name = ('Sample4B(1)-1_CPMG_352015.tnt');
data(26).name = ('Sample4B(1)-2_CPMG_352015.tnt');
data(27).name = ('Sample6B(1)-1_CPMG_352015.tnt');
data(28).name = ('Sample6B(1)-2_CPMG_352015.tnt');
data(29).name = ('Sample8B(1)-1_CPMG_352015.tnt');
data(30).name = ('Sample8B(1)-2_CPMG_352015.tnt');
data(31).name = ('Sample10B(1)-1_CPMG_352015.tnt');
data(32).name = ('Sample10B(1)-2_CPMG_352015.tnt');
data(33).name = ('Sample15B(1)-1_CPMG_352015.tnt');
data(34).name = ('Sample15B(1)-2_CPMG_352015.tnt');
data(35).name = ('Sample20B(1)-1_CPMG_352015.tnt');
data(36).name = ('Sample20B(1)-2_CPMG_352015.tnt');
% the above samples are for the pozzolons post cure


% % these are other samples below, for 0-2% brick dust in 0.25 increments
% data(37).name = ('Sample0_10292014_256scans_512Echoes_CPMG.tnt');
% data(38).name = ('Sample0_10282014_256scans_512Echoes_CPMG.tnt');
% data(39).name = ('Sample1_10202014_256scans_512Echoes_CPMG.tnt');
% data(40).name = ('Sample1_10212014_256scans_512Echoes_CPMG.tnt');
% data(41).name = ('Sample2_10212014_256scans_512Echoes_CPMG.tnt');
% data(42).name = ('Sample2_10222014_256scans_512Echoes_CPMG.tnt');
% data(43).name = ('Sample3_10222014_256scans_512Echoes_CPMG.tnt');
% data(44).name = ('Sample3_10232014_256scans_512Echoes_CPMG.tnt');
% data(45).name = ('Sample4_10232014_256scans_512Echoes_CPMG.tnt');
% data(46).name = ('Sample4_10242014_256scans_512Echoes_CPMG.tnt');
% data(47).name = ('Sample5_10252014_256scans_512Echoes_CPMG.tnt');
% data(48).name = ('Sample5_10242014_256scans_512Echoes_CPMG.tnt');
% data(48).name = ('Sample6_10252014_256scans_512Echoes_CPMG.tnt');
% data(49).name = ('Sample6_10262014_256scans_512Echoes_CPMG.tnt');
% data(50).name = ('Sample7_10262014_256scans_512Echoes_CPMG.tnt');
% data(51).name = ('Sample7_10272014_256scans_512Echoes_CPMG.tnt');
% data(52).name = ('Sample8_10272014_256scans_512Echoes_CPMG.tnt');
% data(53).name = ('Sample8_10282014_256scans_512Echoes_CPMG.tnt');
% there is some complicated indexing below, dependent on the lengths of the
% files, so it would be difficult to add more files to this script. 

% data(37).name = ('Sample1B(1)_CPMG_3192015.tnt');
% data(38).name = ('Sample6B(1)_CPMG_3192015.tnt');

%% reading data in, reshaping the data and cutting out blank points
for i = 1:length(data)
    [data(i).params,data(i).data] = readTecmag4d(data(i).name);
    data(i).data = reshape(real(data(i).data),nrPts,nEchoes);
    data(i).data = data(i).data(1:end-nrBlank,:);
end

% taking an average value from every 3 data files, to use later for the
% fitting
for i=1:length(data)% 
%     sample(i).data= (data(i).data + data(i+1).data); % takes the average
    sample(i).data = data(i).data;% conversion, so averaging doesn't have to take place here.
    sample(i).sumData = sum(sample(i).data);
    sample(i).normSum = sample(i).sumData(:)./max(sample(i).sumData);
end

%% fitting with nlinfit and such
echoAxis = (1:nEchoes)*tE;
guesses = [0.8,0.04,0.2,0.005];
guesses2 = [0.8,0.2];

for i = 1:length(sample)
    [sample(i).beta,sample(i).resid,sample(i).j] = nlinfit(echoAxis(3:end),sample(i).normSum(3:end)',@t2bifit_simple,guesses);
    sample(i).pred = t2bifit_simple(sample(i).beta,echoAxis(3:end));
    sample(i).ci = nlparci(sample(i).beta,sample(i).resid,'jacobian',sample(i).j);
    
    [fixed(i).beta,fixed(i).resid,fixed(i).j] = nlinfit(echoAxis(3:end),sample(i).normSum(3:end)',@t2bifit_monofixed,guesses2);
    fixed(i).pred = t2bifit_monofixed(fixed(i).beta,echoAxis(3:end));
    fixed(i).ci = nlparci(fixed(i).beta,fixed(i).resid,'jacobian',fixed(i).j);

end

%for compiling fits and finding the average
T2 = sample(1).beta(2);
T22 = sample(1).beta(4);
amp = sample(1).beta(1);
amp2 = sample(1).beta(3);
T2eu = sample(1).ci(2,:);
T22eu = sample(1).ci(4,:);
% FT2 = fixed(1).beta(2);
Famp = fixed(1).beta(1);
Fampe = fixed(1).ci(1,:);
Famp2 = fixed(1).beta(2);
Famp2e = fixed(1).ci(2,:);
% FT2eu = fixed(i).ci(2,:);

for i = 2:length(sample)
    T2 = [T2;sample(i).beta(2)];
    T22 = [T22;sample(i).beta(4)];
    amp = [amp; sample(i).beta(1)];
    amp2 = [amp2; sample(i).beta(3)];
    T2eu = [T2eu; sample(i).ci(2,:)];
    T22eu = [T22eu; sample(i).ci(4,:)];
%     FT2 = [FT2;fixed(i).beta(2)];% for fixed values
    Famp = [Famp; fixed(i).beta(1)];
    Fampe = [Fampe; fixed(i).ci(1,:)];
    Famp2 = [Famp2; fixed(i).beta(2)];
    Famp2e = [Famp2e; fixed(i).ci(2,:)];
%     FT2eu = [FT2eu; fixed(i).ci(2,:)];

end

%error stuff to get it in the right format
T2eu = abs([T2eu(:,1) - T2,T2eu(:,2) - T2]);
T22eu = abs([T22eu(:,1) - T22,T22eu(:,2) - T22]);
% FT2eu = abs([FT2eu(:,1) - FT2,FT2eu(:,2) - FT2]);% for fixed fitting 
Fampe = abs([Fampe(:,1) - Famp,Fampe(:,2) - Famp]);
Famp2e = abs([Famp2e(:,1) - Famp2,Famp2e(:,2) - Famp2]);

avgT2 = mean(T2);
avgT22 = mean(T22);
avgAmp = mean(amp);
avgAmp2 = mean(amp2);

%% Averaging the different functions and fits
for i = 1:length(pozzDust)% for pozzolon, this is only indexing to break them apart
    pozzT2sc1(i) = T2(2*i-1);
    pozzT2sc2(i) = T2(2*i);
    pozzT2errorsc1(i,:) = T2eu(2*i-1,:);
    pozzT2errorsc2(i,:) = T2eu(2*i,:);
    
    pozzT22sc1(i) = T22(2*i-1);
    pozzT22sc2(i) = T22(2*i);
    pozzT22errorsc1(i,:,:) = T22eu(2*i-1,:);
    pozzT22errorsc2(i,:) = T22eu(2*i,:);
    
    pozzAmpsc1(i) = amp(2*i-1);
    pozzAmpsc2(i) = amp(2*i);
    pozzAmp2sc1(i) = amp2(2*i-1);
    pozzAmp2sc2(i) = amp2(2*i);
    
    pozzFampsc1(i) = Famp(2*i-1);
    pozzFampsc2(i) = Famp(2*i);
    pozzFamperrorsc1(i,:) = Fampe(2*i-1,:);
    pozzFamperrorsc2(i,:) = Fampe(2*i,:);
    
    pozzFamp2sc1(i) = Famp2(2*i-1);
    pozzFamp2sc2(i) = Famp2(2*i);
    pozzFamp2errorsc1(i,:) = Famp2e(2*i-1,:);
    pozzFamp2errorsc2(i,:) = Famp2e(2*i,:);
end

for i = 1:length(brickDust)% for brickdust
    brickT2sc1(i) = T2(2*i+15);
    brickT2sc2(i) = T2(2*i+16);
    brickT2errorsc1(i,:) = T2eu(2*i+15,:);
    brickT2errorsc2(i,:) = T2eu(2*i+16,:);
    
    brickT22sc1(i) = T22(2*i+15);
    brickT22sc2(i) = T22(2*i+16);
    brickT22errorsc1(i,:) = T22eu(2*i+15,:);
    brickT22errorsc2(i,:) = T22eu(2*i+16,:);
    
    brickAmpsc1(i) = amp(2*i+15);
    brickAmpsc2(i) = amp(2*i+16);
    brickAmp2sc1(i) = amp2(2*i+15);
    brickAmp2sc2(i) = amp2(2*i+16);
    
    brickFampsc1(i) = Famp(2*i+15);
    brickFampsc2(i) = Famp(2*i+16);
    brickFamperrorsc1(i,:) = Fampe(2*i+15,:);
    brickFamperrorsc2(i,:) = Fampe(2*i+16,:);
    
    brickFamp2sc1(i) = Famp2(2*i+15);
    brickFamp2sc2(i) = Famp2(2*i+16);
    brickFamp2errorsc1(i,:) = Famp2e(2*i+15,:);
    brickFamp2errorsc2(i,:) = Famp2e(2*i+16,:);
end
% averaging happens below
pozzAmp = mean([pozzAmpsc1;pozzAmpsc2]);
pozzAmp2 = mean([pozzAmp2sc1;pozzAmp2sc2]);
pozzFamp = mean([pozzFampsc1;pozzFampsc2]);
pozzFamp2 = mean([pozzFamp2sc1;pozzFamp2sc2]);
pozzT2 = mean([pozzT2sc1;pozzT2sc2]);
pozzT22 = mean([pozzT22sc1;pozzT22sc2]);
% and for brick dust

brickAmp = mean([brickAmpsc1;brickAmpsc2]);
brickAmp2 = mean([brickAmp2sc1;brickAmp2sc2]);
brickFamp = mean([brickFampsc1;brickFampsc2]);
brickFamp2 = mean([brickFamp2sc1;brickFamp2sc2]);
brickT2 = mean([brickT2sc1;brickT2sc2]);
brickT22 = mean([brickT22sc1;brickT22sc2]);


%% Plotting all of the functions

figure(1)
plot(pozzDust,pozzAmp,'b','LineStyle','--','LineWidth',1)
hold on
errorbar(pozzDust,pozzFamp,pozzFamperrorsc1(:,1),'b','LineWidth',1.2)% errors are taken from the first scan for convenience, they are very similar though
plot(brickDust,brickAmp,'r','LineStyle','--','LineWidth',1)
errorbar(brickDust,brickFamp,brickFamperrorsc1(:,1),'r','LineWidth',1.2)% errors are taken from the first scan for convenience, they are very similar though

plot(pozzDust,pozzAmp2,'b','LineStyle','--','LineWidth',1)
errorbar(pozzDust,pozzFamp2,pozzFamp2errorsc1(:,1),'b','LineWidth',1.2)% errors are taken from the first scan for convenience, they are very similar though
plot(brickDust,brickAmp2,'r','LineStyle','--','LineWidth',1)
errorbar(brickDust,brickFamp2,brickFamp2errorsc1(:,1),'r','LineWidth',1.2)% errors are taken from the first scan for convenience, they are very similar though

scatter(pozzDust,pozzAmpsc1,40,'b')
scatter(pozzDust,pozzAmpsc2,40,'b')
scatter(pozzDust,pozzAmp2sc1,40,'b')
scatter(pozzDust,pozzAmp2sc2,40,'b')% pozz amplitudes for free biexp fit

scatter(pozzDust,pozzFampsc1,40,'b','+')
scatter(pozzDust,pozzFampsc2,40,'b','+')
scatter(pozzDust,pozzFamp2sc1,40,'b','+')
scatter(pozzDust,pozzFamp2sc2,40,'b','+')%Pozz amp for fixed T2 biexp fit

scatter(brickDust,brickAmpsc1,40,'r')
scatter(brickDust,brickAmpsc2,40,'r')
scatter(brickDust,brickAmp2sc1,40,'r')
scatter(brickDust,brickAmp2sc2,40,'r')

scatter(brickDust,brickFampsc1,40,'r','+')
scatter(brickDust,brickFampsc2,40,'r','+')
scatter(brickDust,brickFamp2sc1,40,'r','+')
scatter(brickDust,brickFamp2sc2,40,'r','+')

xlabel('percentage additive')
ylabel('Normalized amp')
legend('Pozzolon free bifit','Pozzolon fixed T_2 bifit','Brickdust free bifit','Brickdust fixed T_2 bifit')

figure(2)
errorbar(pozzDust,pozzT2,pozzT2errorsc1(:,1),'b','LineWidth',1.2)% errors are taken from the first scan for convenience, they are very similar though
hold on
errorbar(brickDust,brickT2,brickT2errorsc1(:,1),'r','LineWidth',1.2)% errors are taken from the first scan for convenience, they are very similar though
errorbar(pozzDust,pozzT22,pozzT22errorsc1(:,1),'b','LineWidth',1.2)% errors are taken from the first scan for convenience, they are very similar though
errorbar(brickDust,brickT22,brickT22errorsc1(:,1),'r','LineWidth',1.2)% errors are taken from the first scan for convenience, they are very similar though

scatter(pozzDust,pozzT2sc1,40,'b')
scatter(pozzDust,pozzT2sc2,40,'b')
scatter(pozzDust,pozzT22sc1,40,'b')
scatter(pozzDust,pozzT22sc2,40,'b')% pozz amplitudes for free biexp fit

scatter(brickDust,brickT2sc1,40,'r')
scatter(brickDust,brickT2sc2,40,'r')
scatter(brickDust,brickT22sc1,40,'r')
scatter(brickDust,brickT22sc2,40,'r')% pozz amplitudes for free biexp fit

legend('Pozzolon T_2','Brickdust T_2')
xlabel('percentage additive')
ylabel('T_2')

