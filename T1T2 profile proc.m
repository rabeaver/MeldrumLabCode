clear
close all
clc

%% User parameters/Load data

% Scripts folder containing T1/T2 fitting routines
% addpath('/Users/tyler/Documents/Research/Aachen/Data/Testing/Scripts/')
% addpath('/Users/tyler/Dropbox/Coding/Matlab Processing/Scripts/')

cd('/Users/tyler/Dropbox/Data/Paint/Forged_Paintings_Berlin/Beltracci_27Feb2012/Seinebruecke/7/')
parfilestem = 'acqu';

painting_title = 'Musikant';
painting_artist = 'von der Schlichten';
painting_year = 1731;
painting_id = 'music';
biexp_T2fit = 0;
T2rough_fit = 1;
omit_point = 1;
SS_off = 0;

opts = statset('nlinfit');
opts.TolFun = 1e-15;
opts.TolX = 1e-15;

%%
params.acqTime = readpar_Kea(strcat(parfilestem,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(parfilestem,'.par'),'bandwidth');
params.nrScansT1 = readpar_Kea(strcat(parfilestem,'.par'),'nrScans');
params.nrScansT2 = readpar_Kea(strcat(parfilestem,'.par'),'nrScansT2');
params.rxPhase = readpar_Kea(strcat(parfilestem,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(parfilestem,'.par'),'rxGain');
params.nrPtsT1 = readpar_Kea(strcat(parfilestem,'.par'),'nrPntsT1');
params.t1Max = readpar_Kea(strcat(parfilestem,'.par'),'tMax');
params.t1Est = readpar_Kea(strcat(parfilestem,'.par'),'t1Est');
params.repTimeT1 = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.repTimeT2 = readpar_Kea(strcat(parfilestem,'.par'),'repTimeT2');
params.b1Freq = readpar_Kea(strcat(parfilestem,'.par'),'b1Freq');
params.stepSize = readpar_Kea(strcat(parfilestem,'.par'),'stepSize');
params.finalDepth = readpar_Kea(strcat(parfilestem,'.par'),'finalDepth');
params.initDepth = readpar_Kea(strcat(parfilestem,'.par'),'initDepth');
params.nrEchoesT2 = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoesT2');
params.nrEchoes = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoes');

% Load data
try
    depthInfo = dlmread('depthInfo.dat');
    depth = depthInfo(:,2);
    T2catch = depthInfo(:,3);
    loadIf = depthInfo(:,4);
    noRough = 0;
catch EM
    EM
    noRough = 1;
    try
        depth = dlmread('depths.dat');
    catch
        if params.finalDepth == params.initDepth
            depth = params.finalDepth;
        else
           error('proc:BadDepthInfo','The depth info provided by Prospa is incorrect.')
        end
    end
    loadIf = ones(length(depth),1);
end

%% Load T1 and T2 data for detailed measurements
m=1;

if max(loadIf ~= 0)
    
    for n = 1:1:size(loadIf)
        if loadIf(n) == 1
            v = genvarname(['T1data' int2str(depth(n))]);
            w = genvarname(['T2data' int2str(depth(n))]);
            eval([v ' = dlmread(''T1data-' num2str(depth(n)) '.dat'');']);
            eval([w ' = dlmread(''T2data-' num2str(depth(n)) '.dat'');']);
            measDepth(m) = depth(n);
            m = m+1;
        end
    end
    
    allT1Re = zeros(params.nrPtsT1,size(measDepth,2));
    allT1Im = zeros(params.nrPtsT1,size(measDepth,2));
    allT1Mag = zeros(params.nrPtsT1,size(measDepth,2));
    allT2Re = zeros(params.nrEchoesT2,size(measDepth,2));
    allT2Im = zeros(params.nrEchoesT2,size(measDepth,2));
    
    for n = 1:1:size(measDepth,2)
        T1time = eval(['T1data' int2str(measDepth(n)) '(:,1);']);
        T2time = eval(['T2data' int2str(measDepth(n)) '(:,1);']);
        allT1Re(:,n) = eval(['T1data' int2str(measDepth(n)) '(:,2);']);
        allT1Im(:,n) = eval(['T1data' int2str(measDepth(n)) '(:,3);']);
        allT1Mag(:,n) = eval(['sqrt(T1data' int2str(measDepth(n)) '(:,3).^2 + T1data' int2str(measDepth(n)) '(:,2).^2);']);
        allT2Re(:,n) = eval(['T2data' int2str(measDepth(n)) '(:,2);']);
        allT2Im(:,n) = eval(['T2data' int2str(measDepth(n)) '(:,3);']);
    end
    
end

if (noRough ~= 1 && SS_off ~= 1);
    for n = 1:1:size(depth,1)
        w = genvarname(['T2data_rough_' int2str(depth(n))]);
        if ispc == 1
           eval([w ' = dlmread(''rough\T2data_rough-' num2str(depth(n)) '.dat'');']); 
        else
            eval([w ' = dlmread(''./rough/T2data_rough-' num2str(depth(n)) '.dat'');']);
        end
        T2rough(:,n) = eval(['complex(T2data_rough_' int2str(depth(n)) '(:,2),T2data_rough_' int2str(depth(n)) '(:,2));']);
        [T2roughPhased(:,n),~] = autophase(T2rough(:,n),0.1);
    end
elseif (noRough ~= 1 && SS_off == 1);
    for n = 1:1:size(depth,1)
        T2rough(:,n) = eval(['complex(T2data' int2str(depth(n)) '(:,2),T2data' int2str(depth(n)) '(:,2));']);
        [T2roughPhased(:,n),~] = autophase(T2rough(:,n),0.1);
    end
end

if m == 1;
    T2time = eval(['T2data_rough_' int2str(depth(1)) '(:,1);']);
end

%% Plot all data
if m ~= 1;
    figure(1)
    hold all
    h(1) = subplot(2,2,1);
    meshc(measDepth,T1time,allT1Re)
    xlabel('Depth (um)')
    ylabel('time (ms)')
    zlabel('signal intensity (arb)')
    title('T1 Real')
    h(2) = subplot(2,2,2);
    meshc(measDepth,T1time,allT1Im)
    xlabel('Depth (um)')
    ylabel('time (ms)')
    zlabel('signal intensity (arb)')
    title('T1 Imag')
    h(3) = subplot(2,2,3);
    meshc(measDepth,T2time(omit_point+1:end),allT2Re(omit_point+1:end,:))
    xlabel('Depth (um)')
    ylabel('time (ms)')
    zlabel('signal intensity (arb)')
    title('T2 Real')
    h(4) = subplot(2,2,4);
    meshc(measDepth,T1time,allT1Mag)
    xlabel('Depth (um)')
    ylabel('time (ms)')
    zlabel('signal intensity (arb)')
    title('T1 Mag')
end

%% Detailed T2 plot
figure(2)
plot(depth,real(T2roughPhased(omit_point+1,:)))

figure(3)
mesh(depth,T2time(omit_point+1:end),real(T2roughPhased(omit_point+1:end,:)))

%% Fit rough data to T2 decay

if T2rough_fit == 1
   
    y0_guess = 0;
    A_guess = max(max(T2roughPhased));
    t2_guess = T2time(5); %ms
    
    if biexp_T2fit == 1;
        guesses = [y0_guess;A_guess;t2_guess/10;A_guess/2;t2_guess];
        
        CI = 90; %desired confidence interval in percent
        
        ypred_T2_rough = zeros(1001,length(depth));
        
        for n = 1:1:length(depth)
            try
                [xfit_T2,ypred_T2_rough(:,n),T2_coeffs] = bidecay_t2fit(T2time(omit_point+1:end),real(T2roughPhased(omit_point+1:end,n)),guesses,CI);
                y0_T2_rough(n,:) = T2_coeffs(1,1:2);
                A_T2_rough(n,:) = T2_coeffs(2,1:2);
                tau_T2_rough(n,:) = T2_coeffs(3,1:2);
                A2_T2_rough(n,:) = T2_coeffs(4,1:2);
                tau2_T2_rough(n,:) = T2_coeffs(5,1:2);
            catch
                y0_T2_rough(n,:) = [0;0];
                A_T2_rough(n,:) = [0;0];
                tau_T2_rough(n,:) = [0;0];
                A2_T2_rough(n,:) = [0;0];
                tau2_T2_rough(n,:) = [0;0];
            end
        end
        
        tau_T2_rough(:,3) = tau_T2_rough(:,2)./tau_T2_rough(:,1)*100;
        tau2_T2_rough(:,3) = tau2_T2_rough(:,2)./tau2_T2_rough(:,1)*100;
        
    else
        guesses = [y0_guess;A_guess;t2_guess];
        
        CI = 90; %desired confidence interval in percent
        
        ypred_T2_rough = zeros(1001,length(depth));
        
        for n = 1:1:length(depth)
            try
                [xfit_T2,ypred_T2_rough(:,n),T2_coeffs] = monodecay_t2fit(T2time(omit_point+1:end),real(T2roughPhased(omit_point+1:end,n)),guesses,CI,opts);
                y0_T2_rough(n,:) = T2_coeffs(1,1:2);
                A_T2_rough(n,:) = T2_coeffs(2,1:2);
                tau_T2_rough(n,:) = T2_coeffs(3,1:2);
            catch
                y0_T2_rough(n,:) = [0;0];
                A_T2_rough(n,:) = [0;0];
                tau_T2_rough(n,:) = [0;0];
            end
        end
        
        tau_T2_rough(:,3) = tau_T2_rough(:,2)./tau_T2_rough(:,1)*100;
    end
    
end

%% Fit to T2 decay

y0_guess = 0;
A_guess = max(max(allT2Re));
t2_guess = T2time(5); %ms

if biexp_T2fit == 1;
    guesses = [y0_guess;A_guess;t2_guess/10;A_guess/2;t2_guess];
    
    CI = 90; %desired confidence interval in percent
    
    ypred_T2 = zeros(1001,size(measDepth,2));
    
    for n = 1:1:size(measDepth,2)
        [xfit_T2,ypred_T2(:,n),T2_coeffs] = bidecay_t2fit(T2time(omit_point+1:end),allT2Re(omit_point+1:end,n),guesses,CI);
        y0_T2(n,:) = T2_coeffs(1,1:2);
        A_T2(n,:) = T2_coeffs(2,1:2);
        tau_T2(n,:) = T2_coeffs(3,1:2);
        A2_T2(n,:) = T2_coeffs(4,1:2);
        tau2_T2(n,:) = T2_coeffs(5,1:2);
    end
    
    tau_T2(:,3) = tau_T2(:,2)./tau_T2(:,1)*100;
    tau2_T2(:,3) = tau2_T2(:,2)./tau2_T2(:,1)*100;

else 
    guesses = [y0_guess;A_guess;t2_guess]; 
    
    CI = 90; %desired confidence interval in percent
    
    ypred_T2 = zeros(1001,size(measDepth,2));
    
    for n = 1:1:size(measDepth,2)
        [xfit_T2,ypred_T2(:,n),T2_coeffs(:,n),T2_error(:,n),T2_resid(:,n)] = monodecay_t2fit(T2time(omit_point+1:end),allT2Re(omit_point+1:end,n),guesses,CI,opts);
        y0_T2(n) = T2_coeffs(1,n);
        A_T2(n) = T2_coeffs(2,n);
        tau_T2(n) = T2_coeffs(3,n);
    end
    
    tau_T2(:,3) = tau_T2(:,2)./tau_T2(:,1)*100;
end

%% Fit to T1 decay

y0_guess = 0;
A_guess = max(max(allT1Re));
t1_guess = T1time(4); %ms

guesses = [y0_guess;A_guess;t1_guess];

CI = 90; %desired confidence interval in percent

ypred_T1 = zeros(1001,size(measDepth,2));

for n = 1:1:size(measDepth,2)
    [xfit_T1,ypred_T1(:,n),T1_coeffs] = T1fit(T1time,allT1Re(:,n),guesses,CI);
    y0_T1(n,:) = T1_coeffs(1,1:2);
    A_T1(n,:) = T1_coeffs(2,1:2);
    tau_T1(n,:) = T1_coeffs(3,1:2);
end

tau_T1(:,3) = tau_T1(:,2)./tau_T1(:,1)*100;

%% Save data to file

% fid = fopen('C:\Users\TMeldrum\Dropbox\Data\Pinakotheken_May2012\DataForDatabase.txt', 'a');
% for i = 1:1:length(measDepth)
%     fprintf(fid, '%s\t %s\t %4i\t %s\t %6f\t %6f\t %6f\t %6f\t\n', painting_title, painting_artist, painting_year, painting_id, 1000*tau_T1(i,1), 1000*tau_T1(i,2), 1000*tau_T2(i,1), 1000*tau_T2(i,2));
% end
% fclose(fid);

%% Plot fitted data

figure(4)
hold all
h(1) = subplot(2,2,1);
meshc(measDepth,xfit_T1.*1000,ypred_T1)
xlabel('Depth (um)')
ylabel('time (ms)')
zlabel('signal intensity (arb)')
title('T1 Fit')
h(2) = subplot(2,2,2);
meshc(measDepth,xfit_T2,ypred_T2)
xlabel('Depth (um)')
ylabel('time (ms)')
zlabel('signal intensity (arb)')
title('T2 Fit')
h(3) = subplot(2,2,3);
errorbar(measDepth,tau_T1(:,1).*1000,tau_T1(:,2).*1000,'ob')
xlabel('Depth (um)')
ylabel('T1 (ms)')
h(4) = subplot(2,2,4);
hold on
errorbar(measDepth,tau_T2(:,1).*1000,tau_T2(:,2).*1000,'ob')
if biexp_T2fit == 1;
    errorbar(measDepth,tau2_T2(:,1).*1000,tau2_T2(:,2).*1000,'or')
%     legend('T_2 1','T_2 2')
end
xlabel('Depth (um)')
ylabel('T2 (us)')

%% T1 scaling and refitting by group

for n = 1:1:size(measDepth,2)
    t1time_scaled(:,n) = T1time;
    t1data_scaled(:,n) = (allT1Mag(:,n)+y0_T1(n,1))./A_T1(n,1);
end

t1time_scaled = reshape(t1time_scaled,size(t1time_scaled,1)*size(t1time_scaled,2),1);
t1data_scaled = reshape(t1data_scaled,size(t1data_scaled,1)*size(t1data_scaled,2),1);

scaledT1 = [t1time_scaled,t1data_scaled];
scaledT1 = sortrows(scaledT1);

y0_guess = 0;
A_guess = max(scaledT1(:,2));
t1_guess = scaledT1(5*size(measDepth,2)); %ms

guesses = [y0_guess;A_guess;t1_guess];

CI = 90; %desired confidence interval in percent

[xfit_T1scaled,ypred_T1scaled,T1_coeffsscaled] = T1fit(scaledT1(:,1),scaledT1(:,2),guesses,CI);
    y0_T1scaled = T1_coeffsscaled(1,1:2);
    A_T1scaled = T1_coeffsscaled(2,1:2);
    tau_T1scaled = T1_coeffsscaled(3,1:2);


tau_T1scaled(:,3) = tau_T1scaled(:,2)./tau_T1scaled(:,1)*100;

figure
hold on
scatter(scaledT1(:,1),scaledT1(:,2))
plot(xfit_T1scaled,ypred_T1scaled,'-k')

%% T2 scaling and refitting by group

for n = 1:1:size(measDepth,2)
    t2time_scaled(:,n) = T2time;
    t2data_scaled(:,n) = (allT2Re(:,n)+y0_T2(n,1))./A_T2(n,1);
end

t2time_scaled = reshape(t2time_scaled,size(t2time_scaled,1)*size(t2time_scaled,2),1);
t2data_scaled = reshape(t2data_scaled,size(t2data_scaled,1)*size(t2data_scaled,2),1);

scaledT2 = [t2time_scaled,t2data_scaled];
scaledT2 = sortrows(scaledT2);
scaledT2 = scaledT2(5:end,:);

y0_guess = 0;
A_guess = max(scaledT2(:,2));
t2_guess = scaledT2(20*size(measDepth,2))/10; %ms
A_guess2 = max(scaledT2(:,2))/2;
t2_guess2 = scaledT2(20*size(measDepth,2)); %ms

guesses = [y0_guess;A_guess;t2_guess;A_guess2;t2_guess2];

CI = 90; %desired confidence interval in percent

[xfit_T2scaled,ypred_T2scaled,T2_coeffsscaled] = bidecay_t2fit(scaledT2(:,1),scaledT2(:,2),guesses,CI);
    y0_T2scaled = T2_coeffsscaled(1,1:2);
    A_T2scaled = T2_coeffsscaled(2,1:2);
    tau_T2scaled(1,:) = T2_coeffsscaled(3,1:2);
    tau_T2scaled(2,:) = T2_coeffsscaled(5,1:2);


tau_T2scaled(:,3) = tau_T2scaled(:,2)./tau_T2scaled(:,1)*100;

figure
hold on
scatter(scaledT2(:,1),scaledT2(:,2))
plot(xfit_T2scaled,ypred_T2scaled,'-k')

