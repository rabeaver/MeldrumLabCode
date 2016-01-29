clear
clc
close all

%%
data = load('/Users/tyler/Desktop/6Feb2015/T2_D/5/dataRe.dat');

figure(1)
surf(data)
shading flat


gamma = 2.675222005e8; %s-1 T-1
gammaM = 42.57748; %MHz T-1
GM = 1016; %MHz m-1
G = GM/gammaM; %T m-1
DELTA = 10e-3; %s
%%
s = data(:,1);
tau = logspace(log10(15e-6),log10(1015e-6),41);
s1 = log10(s./max(s));
s_x = -gamma^2*G^2.*tau.^2.*(DELTA+tau./3);

plot(s)
%%
% plot(s_x,s1)

plot(tau.^2.*(DELTA+tau/3),s1)


% ln(s/s0) = -gamma^2*D*G^2*tau^2*(DELTA+tau/3)
%            (s^-2 T^-2 m^2 s^-1 T^2 m^-2 s^3)
