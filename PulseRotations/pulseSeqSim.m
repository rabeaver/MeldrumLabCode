% Simluate pulses and precession
% TKM, 5 May 2014

clear
clc
close all

%%
% Initial Magnetization

M(:,1) = [1/sqrt(2);    %Mx
          1/sqrt(2);    %My
          0];   %Mz
  
delw0 = 10;      %Frequency offset

%%
%90 Pulse

%calibration for 90 pulse

flip = 90;                      %flip angle (degrees)
t_cal = 3e-6;                   %pulse length
w1_90 = sqrt((flip*(2*pi/360)/t_cal)^2-delw0^2); %pulse strength (90)
flip = 180;                      %flip angle (degrees)
w1_180 = sqrt((flip*(2*pi/360)/t_cal)^2-delw0^2); %pulse strength (90)

tp=3e-6; %pulse length in s
ph_d = 0; %phase (0 = x, 90 = y, 180 = -x, 270 = -y)

tE = 100e-3;                %interpulse delay
t = 0:1e-5:tE;              %interpulse time axis

M(:,2) = pulseRotation(M(:,1),tp,delw0,ph_d,w1_90); %90
M(:,3) = freePrecession(M(:,2),45e-6,delw0);
M(:,4) = pulseRotation(M(:,3),tp,delw0,ph_d,w1_180);%180
for i = 1:length(t);
% M(:,i) = pulseRotation(M(:,i-1),t,delw0,ph_d,w1);
Mp(:,i) = freePrecession(M(:,4),t(i),delw0);
end
M(:,5) = Mp(:,length(t));

figure(1)
hold on
plot(t,real(Mp(1,:)),'-k')
plot(t,imag(Mp(2,:)),'-r')
plot(t,real(Mp(3,:)),'-b')
ylim([-1 1])
legend('x','y','z')

M(:,6) = pulseRotation(M(:,5),tp,delw0,ph_d,w1_180);

for i = 1:length(t);
% M(:,i) = pulseRotation(M(:,i-1),t,delw0,ph_d,w1);
Mp(:,i) = freePrecession(M(:,6),t(i),delw0);
end


M(:,7) = Mp(:,length(t));

figure(2)
hold on
plot(t,real(Mp(1,:)),'-k')
plot(t,imag(Mp(2,:)),'-r')
plot(t,real(Mp(3,:)),'-b')
ylim([-1 1])
legend('x','y','z')
