close all
clear
clc
%% setting up parameters

nEchoes = 512
tE = 150e-6
nrPts = 69
nrBlank = 5 
pozzDust = [0 0.5 1 2 3 4 5 10];% pozzolon in samples
brickDust = [1 2 3 4 6 8 10 15 20];% brick dust in samples
nr2Dpts = 7;
%% setting up directories

dir1 = ('C:\Users\benjamin\Documents\Data\Water Saturation B-P range\'); % personal computer
cd(dir1)

%% reading in data

data(1).name = ('Sample0P(2)_6.5ml_CPMG_472015.tnt');
data(2).name = ('Sample0.5P(2)_6.5ml_CPMG_472015.tnt');
data(3).name = ('Sample1P(2)_6.5ml_CPMG_472015.tnt');
data(4).name = ('Sample2P(2)_6.5ml_CPMG_4102015.tnt');
data(5).name = ('Sample3P(2)_6.5ml_CPMG_472015.tnt');
data(6).name = ('Sample4P(2)_6.5ml_CPMG_472015.tnt');
data(7).name = ('Sample5P(2)_6.5ml_CPMG_472015.tnt');
data(8).name = ('Sample10P(2)_6.5ml_CPMG_472015.tnt');
% pozzolons above, brick dust below
data(9).name = ('Sample1B(2)_6.5ml_CPMG_482015.tnt');
data(10).name = ('Sample2B(2)_6.5ml_CPMG_492015.tnt');
data(11).name = ('Sample3B(2)_6.5ml_CPMG_492015.tnt');
data(12).name = ('Sample4B(2)_6.5ml_CPMG_492015.tnt');
data(13).name = ('Sample6B(2)_6.5ml_CPMG_492015.tnt');
data(14).name = ('Sample8B(2)_6.5ml_CPMG_492015.tnt');
data(15).name = ('Sample10B(2)_6.5ml_CPMG_492015.tnt');
data(16).name = ('Sample15B(2)_6.5ml_CPMG_492015.tnt');
data(17).name = ('Sample20B(2)_6.5ml_CPMG_4102015.tnt');

%% reading data in, reshaping the data and cutting out blank points
for i = 1:length(data)
    [data(i).params,~,data(i).twoD] = readTecmag4d(data(i).name); % reads it into parameters
end

for j =1:(length(data)) % removes the 2d Matrix structures for ease of processing
    for i = 1:nr2Dpts
    sample(i+(j-1)*nr2Dpts).data = data(j).twoD(i,:);
    end
end

for i = 1:length(sample)% reshapes
    sample(i).data = reshape(real(sample(i).data),nrPts,nEchoes);
    sample(i).data = sample(i).data(1:end-nrBlank,:);% and cuts out blank points
    sample(i).sumData = sum(sample(i).data);
    sample(i).normSum = sample(i).sumData(:)./max(sample(i).sumData); % normalizes data
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

T2points = reshape(T2,nr2Dpts,length(data));
T22points = reshape(T22,nr2Dpts,length(data));
ampPoints = reshape(amp,nr2Dpts,length(data));
amp2Points = reshape(amp2,nr2Dpts,length(data));
Famp = reshape(Famp,nr2Dpts,length(data));
Famp2 = reshape(Famp2,nr2Dpts,length(data));

T2euPoints = reshape(T2eu,nr2Dpts,length(data),2);
T22euPoints = reshape(T22eu,nr2Dpts,length(data),2);% creates a 3d matrix, w/(~,~,1) being the 1st and (~,~,2) being the 2nd error point for each corresponding point
Fampe = reshape(Fampe,nr2Dpts,length(data),2);
Famp2e = reshape(Famp2e,nr2Dpts,length(data),2);

MeanT2points = mean(T2points);
MeanT22points = mean(T22points);
MeanAmpPoints = mean(ampPoints);
MeanAmp2Points = mean(amp2Points);
MeanFamp = mean(Famp);
MeanFamp2 = mean(Famp2);


%% Plotting the functions and such

figure(1)
plot(pozzDust,MeanT2points(1:length(pozzDust)),'b')
hold on
plot(brickDust,MeanT2points((length(pozzDust)+1):end),'r')
plot(pozzDust,MeanT22points(1:length(pozzDust)),'b')
plot(brickDust,MeanT22points((length(pozzDust)+1):end),'r')
legend('Pozzolon','Brickdust')
for i = 1:length(pozzDust)
    for j = 1:nr2Dpts
    scatter(pozzDust(i),T2points(j,i),'b')
    scatter(pozzDust(i),T22points(j,i),'b')
    end
end
for i = 1:length(brickDust)
    for j = 1:nr2Dpts
    scatter(brickDust(i),T2points(j,i+length(pozzDust)),'r')
    scatter(brickDust(i),T22points(j,i+length(pozzDust)),'r')
    end
end
title('T_2 plots')
xlabel('Percentage additive')
ylabel('T_2 [s]')

figure(2)
plot(pozzDust,MeanAmpPoints(1:length(pozzDust)),'b')
hold on
plot(brickDust,MeanAmpPoints((length(pozzDust)+1):end),'r')
% for the fixed amplitude fit below
plot(pozzDust,MeanFamp2(1:length(pozzDust)),'g','LineStyle','--')
plot(brickDust,MeanFamp2((length(pozzDust)+1):end),'k','LineStyle','--')
plot(pozzDust,MeanFamp(1:length(pozzDust)),'g','LineStyle','--')
plot(brickDust,MeanFamp((length(pozzDust)+1):end),'k','LineStyle','--')
plot(pozzDust,MeanAmp2Points(1:length(pozzDust)),'b')
plot(brickDust,MeanAmp2Points((length(pozzDust)+1):end),'r')

legend('Pozzolon','Brickdust','Fixed amp Pozzolon','Fixed amp Brickdust')
for i = 1:length(pozzDust)
    for j = 1:nr2Dpts
    scatter(pozzDust(i),ampPoints(j,i),'b')
    scatter(pozzDust(i),amp2Points(j,i),'b')
    scatter(pozzDust(i),Famp(j,i),'g')
    scatter(pozzDust(i),Famp2(j,i),'g')
    end
end
for i = 1:length(brickDust)
    for j = 1:nr2Dpts
    scatter(brickDust(i),ampPoints(j,i+length(pozzDust)),'r')
    scatter(brickDust(i),amp2Points(j,i+length(pozzDust)),'r')
    scatter(brickDust(i),Famp(j,i+length(pozzDust)),'k')
    scatter(brickDust(i),Famp2(j,i+length(pozzDust)),'k')
    end
end
title('Amplitude plot')
xlabel('Percentage additive')
ylabel('Relative amplitude')


