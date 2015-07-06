clear
close all
clc

% before loading data, need to make sure it looks okay in the Matlab text
% editor. Unicode (UTF-8), Unix LF, and get rid of funck characters at the
% beginning of the file

data = load('~/Dropbox/Data/DSC/PlayDSC.dat');

time_min = data(:,1);
temp_C = data(:,2);
heatFlow_mW = data(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%enter in parameter values for the measurement%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mass = 1.58e-3; %g
MM = 66463; %g mol-1
dT_dt = 0.5; %K min-1
alpha = 0.5; %between 0.25 and 1, for baseline fitting, 1 is steeper
dT_point = 0.01;
nRes = 1166;



heatFlow_Wperg = heatFlow_mW/1000/mass;

temp_C_new = floor(min(temp_C)):dT_point:ceil(max(temp_C));
heatFlow_new = interp1(temp_C,heatFlow_Wperg,temp_C_new);


figure(1)
plot(temp_C_new,heatFlow_new)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at the plot, pick a low and high T to integrate over%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_min = 63.3; %C
T_max = 93; %C


%calculate necessary indices for baseline markers
[~,i_trans] = max(heatFlow_new);
T_trans = temp_C_new(i_trans);

i_low = int16((1/dT_point)*(T_min - floor(min(temp_C)))+1);
i_hi  = int16((1/dT_point)*(T_max - floor(min(temp_C)))+1);
heat_low = heatFlow_new(i_low);
heat_hi  = heatFlow_new(i_hi);

% adjust figure limits to match temperature range
figure(1)
xlim([T_min T_max])

%generare, plot, and subtract sigmoidal horizontal baseline
baseline = heat_low + (heat_hi - heat_low)./(1+exp((-alpha)*(temp_C_new-T_trans)));

figure(1)
hold on
plot(temp_C_new,baseline,'-r')
title('heat flow (W g^{-1} and baseline vs. T')
hold off

heat_corr = heatFlow_new - baseline;
Cp_g = heat_corr/(dT_dt/60)/1000; %kJ g-1 K-1

%plot heat capacity per gram
figure(2)
plot(temp_C_new,Cp_g)
xlim([T_min T_max])
title('C_p vs. T')
ylabel('kJ g-1 K-1')

%determine and plot molar heat capacity
Cp_mol = Cp_g*MM;

figure(3)
plot(temp_C_new,Cp_mol)
xlim([T_min T_max])
title('C_{p,mol} vs. T')
ylabel('kJ mol-1 K-1')

%calc Cp/T dT and integrate to get S
%integrate Cp dT to get H

delS = cumsum(Cp_mol(i_low:i_hi)./temp_C_new(i_low:i_hi))*dT_point;
delH = cumsum(Cp_mol(i_low:i_hi))*dT_point;

figure(4)
hold on
plot(temp_C_new(i_low:i_hi),delH)
plot(temp_C_new(i_low:i_hi),delS)
xlim([T_min T_max])


Y = delH-min(delH);
Y = Y./max(Y);

dY_dT = diff(Y)/dT_point;
dY_dT_Tm = dY_dT(i_trans-i_low);

dH_VH = 4*8.3145e-3*(T_trans+273.15)^2*dY_dT_Tm;
dH = max(delH); %kJ mol-1
dS = max(delS); %kJ mol-1 K-1



