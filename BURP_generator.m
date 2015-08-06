% BURP pulse generator
% TKM, 6 August 2015

clear
close all
clc

%%
pulseLength = 50e-6; %s
flipAngle = 0.5; %give as fraction of 2pi. So a 90 deg flip is 0.25, 180 is 0.5, etc.

%%

nPts = 256;
pulseTimeRes = pulseLength/nPts; %s
pulseTime = pulseTimeRes:pulseTimeRes:pulseLength;

% coefficients for iBURP-1
A = [0.70, -0.15, -0.94, 0.11, -0.02, -0.04, 0.01, -0.02, -0.01];
B = [-1.54, 1.01, -0.24, -0.04, 0.08, -0.04, -0.01, 0.01, -0.01];

for i = 1:length(A)
    p(:,i) = A(i)*cos(i*2*pi/pulseLength*pulseTime) + B(i)*sin(i*2*pi/pulseLength*pulseTime);
end
%%
p = sum(p,2) + flipAngle;

plot(pulseTime,p)

