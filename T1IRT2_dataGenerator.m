% 7 July 2015
% This file generates 2D T1(inversion recovery)-T2 data. The user should
% input T1 and T2 values (only one peak is in this implementation), the
% time increment for both the T1 and T2 axes, and the number of points on
% each axis. The program generates T2 and T1 curves, then multiplies them
% (line 32) to make the final data set. Plotting and saving for use in
% Prospa is also avaiable.


clear
close all
clc

%%

% User defined parameters
T1 = 4.5e-3; %s
T2 = 2.0e-3; %s
t_T1 = 500e-6; %s, time increment for T1 axis
t_T2 = 500e-6; %s, time increment for T2 axis
nPtsT1 = 50; %points for T1 axis
nPtsT2 = 200; %points for T2 axis
%end user-defined parameters

data = zeros(nPtsT1,nPtsT2); %make blank matrix
T1axis = logspace(log10(t_T1),log10(t_T1*nPtsT1),nPtsT1); %set up axes
T2axis = logspace(log10(t_T2),log10(t_T2*nPtsT2),nPtsT2);

T2data = exp(-T2axis/T2); %1D T2 curve
T1data = 1-2*exp(-T1axis/T1); %1D T1 curve

data2D = T1data'*T2data; %make 2D data set

surf(T2axis',T1axis,data2D) %plotting
shading flat
set(gca,'Xscale','log','Yscale','log')
ylabel('T1 [s]')
xlabel('T2 [s]')
ylim([min(T1axis) max(T1axis)]);
xlim([min(T2axis) max(T2axis)]);