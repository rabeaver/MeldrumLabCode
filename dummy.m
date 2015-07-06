clc
close all
clear

%% Load all data

datadir = '/Users/tyler/Documents/Research/Aachen/Data/T1 vs age data/BB Paintings/BB_lavallchateau_05mar2012/lowerleft/1/';
autophase = 1;
omitfirstpoint = 1;
T1mag = 1;

% If this is a set of T1/T2 data that has a depthInfo file, it will load
% here:
% [params,allDepth,measDepth,T1time,T2time,T1data,T2data,T2roughdata] = loadFromDepthInfo(datadir,omitfirstpoint,autophase);
[params,measDepth,T1time,T2time,T1data,T2data] = loadFromDepths(datadir,omitfirstpoint,autophase);

usedatapoints = 1:1:length(measDepth);

mesh(measDepth,T1time,abs(T1data))

%% T2 monoexponential fit--fine data
usedatapoints = 1:1:24;

CI = 90;

guesses = [0;real(max(max(T2data)));T2time(4)];
clear T2fitcoeffs;

for i =1:1:length(usedatapoints)
    [T2xfit,T2yfit(:,usedatapoints(i)),T2beta] = monodecay_t2fit(T2time,real(T2data(:,usedatapoints(i))),guesses,CI);
    T2fitcoeffs.tau(usedatapoints(i),:) = T2beta(3,:);
    T2fitcoeffs.A(usedatapoints(i),:) = T2beta(2,:);
    T2fitcoeffs.y0(usedatapoints(i),:) = T2beta(1,:);
end
T2fitcoeffs.tau(:,3) = 100*T2fitcoeffs.tau(:,2)./T2fitcoeffs.tau(:,1);
T2fitcoeffs.A(:,3) = 100*T2fitcoeffs.A(:,2)./T2fitcoeffs.A(:,1);
T2fitcoeffs.y0(:,3) = 100*T2fitcoeffs.y0(:,2)./T2fitcoeffs.y0(:,1);

%% T1 fit

if T1mag == 1
    T1data = abs(T1data);
end

guesses = [0;real(max(max(T1data)));max(T1time)/5];
clear T1fitcoeffs;


for i =1:1:length(usedatapoints)
    [T1xfit,T1yfit(:,usedatapoints(i)),T1beta] = T1_fit(T1time,real(T1data(:,usedatapoints(i))),guesses,CI);
    T1fitcoeffs.tau(usedatapoints(i),:) = T1beta(3,:);
    T1fitcoeffs.A(usedatapoints(i),:) = T1beta(2,:);
    T1fitcoeffs.y0(usedatapoints(i),:) = T1beta(1,:);
end
T1fitcoeffs.tau(:,3) = 100*T1fitcoeffs.tau(:,2)./T1fitcoeffs.tau(:,1);
T1fitcoeffs.A(:,3) = 100*T1fitcoeffs.A(:,2)./T1fitcoeffs.A(:,1);
T1fitcoeffs.y0(:,3) = 100*T1fitcoeffs.y0(:,2)./T1fitcoeffs.y0(:,1);

%% Plot
figure(1)
subplot(3,3,1)
hold all
for i = 1:1:length(usedatapoints)
    scatter(T1time,real(T1data(:,usedatapoints(i))))
    plot(T1xfit,T1yfit(:,usedatapoints(i)))
end
subplot(3,3,2)
hold all
for i = 1:1:length(usedatapoints)
    errorbar(measDepth(usedatapoints(i)),T1fitcoeffs.A(usedatapoints(i),1),T1fitcoeffs.A(usedatapoints(i),2),'o')
end
subplot(3,3,3)
hold all
for i = 1:1:length(usedatapoints)
    errorbar(measDepth(usedatapoints(i)),T1fitcoeffs.tau(usedatapoints(i),1),T1fitcoeffs.tau(usedatapoints(i),2),'o')
end
subplot(3,3,4)
hold all
for i = 1:1:length(usedatapoints)
    scatter(T2time,real(T2data(:,usedatapoints(i))))
    plot(T2xfit,T2yfit(:,usedatapoints(i)))
end
subplot(3,3,5)
hold all
for i = 1:1:length(usedatapoints)
    errorbar(measDepth(usedatapoints(i)),T2fitcoeffs.A(usedatapoints(i),1),T2fitcoeffs.A(usedatapoints(i),2),'o')
end
subplot(3,3,6)
hold all
for i = 1:1:length(usedatapoints)
    errorbar(measDepth(usedatapoints(i)),T2fitcoeffs.tau(usedatapoints(i),1),T2fitcoeffs.tau(usedatapoints(i),2),'o')
end
% subplot(3,3,8:9)
% mesh(allDepth,T2time,real(T2roughdata))

%% Write to file

% finalData = [depth,T2_amp,T2_amp_error,T2_tau,T2_tau_error,T1_amp,T1_amp_error,T1_tau,T1_tau_error];

finalData = zeros(length(usedatapoints),9);

for i = 1:1:length(usedatapoints)
    finalData(i,1) = measDepth(usedatapoints(i));
    finalData(i,2) = T2fitcoeffs.A(usedatapoints(i),1);
    finalData(i,3) = T2fitcoeffs.A(usedatapoints(i),2);
    finalData(i,4) = T2fitcoeffs.tau(usedatapoints(i),1);
    finalData(i,5) = T2fitcoeffs.tau(usedatapoints(i),2);
    finalData(i,6) = T1fitcoeffs.A(usedatapoints(i),1);
    finalData(i,7) = T1fitcoeffs.A(usedatapoints(i),2);
    finalData(i,8) = T1fitcoeffs.tau(usedatapoints(i),1);
    finalData(i,9) = T1fitcoeffs.tau(usedatapoints(i),2);
end

dlmwrite(strcat(datadir,'finalData.txt'),finalData,'\t');