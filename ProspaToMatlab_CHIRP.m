clear
close all
clc

%% T1T2

% Load data
cd('C:\users\jnking01\desktop\'); % cd into desired directory
datafile = 'CHIRP_15mMGdH2O_40mspw_sliceheight350um_Td8u_76pts_512scans_50nsWave_10dB_3Aug2015'; % Datafile
invDat = load(strcat(datafile, '.out')); % Load Prospa File

% Match axes to Prospa Inversion
T2lim = [-4 -1]; % T2 limits
T1lim = [-4 -0]; % T1 limits

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
cd('Z:\JNK\CHIRP\Best_Data_Sets\Ville_Update\') % cd into desired directory
saveas(f, datafile, 'tif') % Save as .tif