% Stimulated Echo Diffusion Delay Calculator
% Use this to find d1, d3 values for the traditional T2-D experiments
% ALL TIMES ARE IN MICROSECONDS FOR TNMR
% TKM, 17 Feb 2017

clear
clc
close all

np = 64; % the number of indirect points
DELTA = 3000; %DELTA diffusion time
tE = 250; %echo time
acqt = 164; %acquisition time
pw = 6; %pulse length
dmin = 12; %minimum delta diffusion time
dmax = 3500; %maximum delta diffusion time

lin = 0; %linearly spaced if 1, logarithmically spaced if 0

%%
if lin==1;
    delta = linspace(dmin,dmax,np);
else
    delta = logspace(log10(dmin),log10(dmax),np);
end

d1 = delta'-pw;
d2 = DELTA-pw;
d3 = d1+tE/2;
d4 = (tE - acqt - pw)/2;

