% SINC amplitude table generator for Tecmag

clear
clc
close all

%%%%%%% User-defined parameters %%%%%%%

slice = 0.080; %mm, slice thickness
G = 6.59; %T m-1, B0 field gradient
tau = 450e-6; %s; pulse length
Nlobes = 0; %# of sinc lobes
dt = 0.5e-6; %s, time per pulse point
offset = 0; %Hz
delay = 2e-6; %s between sinc pulses

%%%%%%% END User-defined parameters %%%%%%%

N = round(tau/dt); %number of points per pulse waveform
gamma = 42.576; %MHz T-1
SW = slice*G*gamma*1000; %Hz
delay = zeros(delay/dt,1);

pulseAmp = -N*dt/2+dt/2:dt:N*dt/2-dt/2;
carrierSig(1,:) = cos(2*pi*pulseAmp*offset(1));
% carrierSig(2,:) = cos(2*pi*pulseAmp*offset(2));
% carrierSig(3,:) = cos(2*pi*pulseAmp*offset(3));

B1(:,1) = sin(SW*pi*pulseAmp)./(SW*pi*pulseAmp).*carrierSig(1,:);
% B1(:,2) = sin(SW*pi*pulseAmp)./(SW*pi*pulseAmp).*carrierSig(2,:);
% B1(:,3) = sin(SW*pi*pulseAmp)./(SW*pi*pulseAmp).*carrierSig(3,:);

B1seq = [B1(:,1)]; % delay; B1(:,2); delay; B1(:,3)];

% alp = 0.43;
% Hamm = (1-alp)+alp*cos(SW*pi*pulseAmp/Nlobes);
% 
% B1sinc = B1.*Hamm;
% 
dlmwrite('SINCout.dat',B1seq,'delimiter',' ');
totalTime = dt*length(B1seq)*1e6

figure(1)
hold on
% plot(pulseAmp,B1sinc,'-k')
plot(B1seq,'-r')
% plot(pulseAmp,carrierSig,':b')