clear
close all
clc

%% T1T2

% Load data
cd('/Users/jaredking/Documents/Classes/Chemistry/Research/Summer2016/TradSSET2_M212Samples_23_28June2016/TradSSET2_M212_Fresh_Overnight_24June2016/2/ILTs/'); % cd into desired directory
datafile = 'TradSSET2_M212_Fresh_Overnight_24June2016'; % Datafile
invDat = load(strcat(datafile, '.out')); % Load Prospa File

% Match axes to Prospa Inversion
T2lim = [-5 -0]; % T2 limits
T1lim = [-12 -7]; % T1 limits

% Make Axes
T2axis = logspace(T2lim(1),T2lim(2), size(invDat,1));
T1axis = logspace(T1lim(1),T1lim(2), size(invDat,2));

% Intensity plot matching the Prospa Figure
f = figure(1); 
hold on
surf(T2axis,T1axis,invDat)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('T2 (s)')
ylabel('T1 (s)')
shading flat
view(2) % Sets correct view for saving

% Save figure to bring into Adobe Illustrator
% cd('Z:\JNK\CHIRP\Best_Data_Sets\Ville_Update\') % cd into desired directory
saveas(f, datafile, 'tif') % Save as .tif