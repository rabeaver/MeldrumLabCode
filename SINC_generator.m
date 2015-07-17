% SINC amplitude table generator for Tecmag

clear
clc
close all

%%%%%%% User-defined parameters %%%%%%%

slice = 0.080; %mm, slice thickness
G = 6.59; %T m-1, B0 field gradient
tau = 450e-6; %s; pulse length
Nlobes = 0; %# of sinc lobes
dt = 1e-6; %s, time per pulse point
offset = 0; %Hz
amp = 3; %db

%%%%%%% END User-defined parameters %%%%%%%

N = round(tau/dt); %number of points per pulse waveform
gamma = 42.576; %MHz T-1
SW = slice*G*gamma*1000; %Hz

pulseAmp = -N*dt/2+dt/2:dt:N*dt/2-dt/2;
carrierSig(1,:) = cos(2*pi*pulseAmp*offset(1));

B1(:,1) = sqrt((sin(SW*pi*pulseAmp)./(SW*pi*pulseAmp).*carrierSig(1,:)).^2)*amp;

dlmwrite('SINC_Amp.dat',B1);
totalTime = dt*length(B1)*1e6;
taxis = linspace(dt, N*dt, N);

figure(1)
hold on
% plot(pulseAmp,B1sinc,'-k')
plot(taxis, B1,'-r')
% plot(pulseAmp,carrierSig,':b')