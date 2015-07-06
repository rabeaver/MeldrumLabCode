% Make synthetic data based on T2-D parameters
% 9 Feb 2015, TKM

clear
close all
clc

%%% User input values here

nEchoes = 128;
tEcho = 60e-6; %s
DELTA = 10e-3; %s
T2 = 10e-3; %s
D = 3e-11; %m^2 s^-1
deltaMin = 15e-6; %s
deltaMax = 1015e-6; %s
n2DPts = 41;
GM = 1016; %MHz m-1 (gradient)

%%% calculated values
gamma = 2.675222005e8; %s-1 T-1
gammaM = 42.57748; %MHz T-1
G = GM/gammaM; %T m-1

%%% make echo and delta vectors, multiply

echoVector = tEcho:tEcho:tEcho*nEchoes;
deltaVector = logspace(log10(deltaMin),log10(deltaMax),n2DPts);

s = exp(-gamma^2*G^2*D.*deltaVector.^2.*(DELTA+deltaVector./3));
echoAttenuation = exp(-echoVector./T2);

data = s'*echoAttenuation;

%%% plot results in surface plot
figure(1)
surf(echoVector,deltaVector,data)
shading flat
xlabel('echo times (s)')
ylabel('delta times (s)')
zlabel('signal intensity (arb)')





