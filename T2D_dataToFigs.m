clear
clc
close all

%give file dir, file name (the *.out file from Prospa export2d), T2 and D limits (should be the same for Naproxed
%stuff), and number of points in inverted data.
datadir = '/Users/tyler/Dropbox/Data/Biosensors/Figures/Normalized Images/Data/';
datalist = {'0_10_out';
            '1_9_out';
            '3_7_out';
            '4_6_out';
            '5_5_out';
            '6_4_out';
            '7_3_out';
            '9_1_out';
            '10_0_out'};

T2lims = [1e-4 1e1];
Dlims = [1e-11 1e-9];
ampNorm = 4e-6; %this is a constant to bring the values of all the intensities closer to 1
        
for sampleID = 1:9;        
datafile = (datalist{sampleID});


%load the data and remove 0 values (replace with NaN)
data = load(strcat(datadir,datafile,'.out'));
% data(data == 0) = NaN;
data = data/ampNorm;
nPts = size(data,1);

%this calculates the axes based on the limits above
Daxis = logspace(log10(Dlims(1)),log10(Dlims(2)),nPts);
T2axis = logspace(log10(T2lims(1)),log10(T2lims(2)),nPts);


%plot the T2D data
figure(sampleID)
set(gcf,'Position',[680 60 800 800])
pcolor(T2axis,Daxis,data)
shading flat
set(gca,'XScale','log','YScale','log')
xlim(T2lims)
ylim(Dlims)
caxis([0 1])
colormap parula
ylabel('\itD\rm [m^2 s^{-1}]')
xlabel('\itT\rm_2 [s]')
% view([0,90])
if sampleID > 1
    axis off
end
if(sampleID == 9)
    colorbar
    set(gcf,'Position',[680 60 880 800])
end
% print(strcat(datalist{sampleID},'_FIG.npg'),'-dpng','-opengl','-r600')
end

%%
clear
clc
close all

%give file dir, file name (the *.out file from Prospa export2d), T2 and D limits (should be the same for Naproxed
%stuff), and number of points in inverted data.
datadir = '/Users/tyler/Dropbox/Data/Biosensors/Figures/Normalized Images/Data/';
datalist = {'0_10_out';
            '1_9_out';
            '3_7_out';
            '4_6_out';
            '5_5_out';
            '6_4_out';
            '7_3_out';
            '9_1_out';
            '10_0_out'};

T2lims = [1e-4 1e1];
Dlims = [1e-11 1e-9];
sampleAxis = [0, 0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1];

ampNorm = 4e-6; %this is a constant to bring the values of all the intensities closer to 1
data = zeros(750,750,9);

for sampleID = 1:9;        
datafile = (datalist{sampleID});


%load the data and remove 0 values (replace with NaN)
datatmp = load(strcat(datadir,datafile,'.out'));
% data(data == 0) = NaN;
data(:,:,sampleID) = datatmp/ampNorm;
nPts = size(data,1);

%this calculates the axes based on the limits above
Daxis = logspace(log10(Dlims(1)),log10(Dlims(2)),nPts);
T2axis = logspace(log10(T2lims(1)),log10(T2lims(2)),nPts);
end

cm = colormap(parula);
% cm(1,:) = [0.8 0.8 0.95];
% cm(2,:) = [0.95 0.95 0.95];
% cm(3,:) = [0.95 0.95 0.95];
% cm(4,:) = [0.95 0.95 0.95];

figure(10);
set(gcf,'Position', [1000    10   1000   950])
slice(T2axis,Daxis,sampleAxis,data,[], [], [0.1,0.3,0.5,0.7,0.9]);
shading flat; 
% alpha(1)
set(gca,'XScale','log','YScale','log','XTick', [1e-4; 1e-3; 1e-2; 1e-1; 1e0],'FontSize',18)
xlim(T2lims)
ylim(Dlims)
zlim([0 1])
% caxis([0 1])
colormap(cm)
ylabel('\itD\rm [m^2 s^{-1}]')
xlabel('\itT\rm_2 [s]')
zlabel('BSA:Naproxen ratio (^1H basis)')
view([23 9.5])
print('ImageCube.tif','-dtiff','-r600')
