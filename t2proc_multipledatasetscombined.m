clear
close all
clc

%% User parameters

% Input data location and filestem, number of profiles, depths per profile,
% and echoes per depth
addpath('/Users/tyler/Documents/Research/Aachen/Data/Museum Ludwig, Jan 2012/Scripts')
datafilestem = 'data';
acqufilestem = 'acqu';

textinfo(1) = {'K. van Dongen, Bildnis Ana'};
textinfo(2) = {'1906-1911'};
textinfo(3) = {'Three CPMG measurements'};
textinfo(4) = {'Position 1, red, first point omitted'};

% textinfo(1) = {'M. Pechstein, Das gruene Sofa'};
% textinfo(2) = {'1910'};
% textinfo(3) = {'Three CPMG measurement(s)'};
% textinfo(4) = {'Position 2, dark green, first point omitted'};

% textinfo(1) = {'M. Pechstein, Das gruene Haus'};
% textinfo(2) = {'1909'};
% textinfo(3) = {'Three CPMG measurement(s)'};
% textinfo(4) = {'Position 6, light green, first point omitted'};

% textinfo(1) = {'M. Pechstein, Zwei Frauenakte im Zimmer'};
% textinfo(2) = {'1909'};
% textinfo(3) = {'Three CPMG measurements combined'};
% textinfo(4) = {'Position 1, Pink, first point omitted'};

% textinfo(1) = {'A. Derain, Blick auf St. Paul-de-Vence'};
% textinfo(2) = {'1910'};
% textinfo(3) = {'Three CPMG measurements combined'};
% textinfo(4) = {'Position 1, dark blue, first point omitted'};

omit_points = 1;


dirpath(1,:) = '/Users/tyler/Documents/Research/Aachen/Data/Museum Ludwig, Jan 2012/17.01.2012/vandongen_ana_pos1_t2/1/';
dirpath(2,:) = '/Users/tyler/Documents/Research/Aachen/Data/Museum Ludwig, Jan 2012/17.01.2012/vandongen_ana_pos1_t2/2/';
dirpath(3,:) = '/Users/tyler/Documents/Research/Aachen/Data/Museum Ludwig, Jan 2012/17.01.2012/vandongen_ana_pos1_t2/3/';
% dirpath(4,:) = '/Users/tyler/Documents/Research/Aachen/Data/Museum Ludwig, Jan 2012/17.01.2012/pechstein_haus_pos5_t2/4/';
% dirpath(1,:) = '/Users/tyler/Documents/Research/Aachen/Data/Museum Ludwig, Jan 2012/17.01.2012/pechstein_haus_pos5_t2/5/';

numDataSets = size(dirpath,1);

%% Data set 

for i = 1:1:numDataSets
    cd(dirpath(i,:))
    
    params.numberProfiles = readpar_Kea(strcat(acqufilestem,'.par'),'nrExp');
    params.finalDepth = readpar_Kea(strcat(acqufilestem,'.par'),'finalDepth');
    params.initDepth = readpar_Kea(strcat(acqufilestem,'.par'),'initDepth');
    params.stepSize = readpar_Kea(strcat(acqufilestem,'.par'),'stepSize');
    params.numberDepths = 1 + (params.finalDepth - params.initDepth)/params.stepSize;
    params.numberEchoes = readpar_Kea(strcat(acqufilestem,'.par'),'nrEchoes');
    params.rxPhase = readpar_Kea(strcat(acqufilestem,'.par'),'rxPhase');
    params.bandwidth = readpar_Kea(strcat(acqufilestem,'.par'),'bandwidth');
    params.pulseLength = readpar_Kea(strcat(acqufilestem,'.par'),'pulseLength');
    params.b1Freq = readpar_Kea(strcat(acqufilestem,'.par'),'b1Freq');
    params.numberScans = readpar_Kea(strcat(acqufilestem,'.par'),'nrScans');
    params.repTime = readpar_Kea(strcat(acqufilestem,'.par'),'repTime');
    params.echoTime = readpar_Kea(strcat(acqufilestem,'.par'),'echoTime0');
    
    % Load data
    data = dlmread(strcat(datafilestem,'.dat'));
    echotime(:,i) = data(omit_points+1:end,1);
    amplitude(:,i) = data(omit_points+1:end,2)./data(omit_points+1,2);    
    
end

combineddata(:,1) = reshape(echotime,size(echotime,1)*size(echotime,2),1);
combineddata(:,2) = reshape(amplitude,size(amplitude,1)*size(amplitude,2),1);
combineddata = sortrows(combineddata);
combineddata = combineddata(numDataSets:end,:);


%% Plot

figure(1)
scatter(combineddata(:,1),combineddata(:,2),'or')

%% Fit to exp decay



y0_guess = 1;
A_guess = 1;
t2_guess = 0.2; %ms

guesses = [y0_guess;A_guess;t2_guess];

CI = 90; %desired confidence interval in percent

% for i = 1:1:size(amplitude,2)
    [xfit,ypred,coeffs] = monodecay_t2fit(combineddata(:,1),combineddata(:,2),guesses,CI);
% end

fitinfo = sprintf('tau = %3.3f +/- %3.3f us (%2.2f%%)',1000*coeffs(3,1),1000*coeffs(3,2),100*coeffs(3,2)/coeffs(3,1));
textinfo(5) = {fitinfo};

figure(2)
hold all
scatter(combineddata(:,1),combineddata(:,2),'ok')
plot(xfit,ypred,'-k')
xlabel('time [ms]')
ylabel('normalized signal amplitude')
text(0.4*max(xfit),0.8*max(ypred),textinfo)

%% Save figure, be sure to get correct params
cd('/Users/tyler/Documents/Research/Aachen/Data/Museum Ludwig, Jan 2012/Final data/');
print -f2 -r150 -depsc vandongen_ana_pos1_combined_t2
