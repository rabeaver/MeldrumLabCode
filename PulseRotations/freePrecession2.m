clear
clc
close all
% Free precession
delw0 = 0.1;
t = 0:0.1:100;


M_0 = [1;   %x
       0;   %y
       0];  %z

for i=1:length(t)
R(:,:,i) = [exp(1i*delw0*t(i))     0            0;
             0              exp(-1i*delw0*t(i)) 0;
             0                  0               1];
 
M_t(:,i) = R(:,:,i)*M_0;
end

figure
hold on
plot(t,real(M_t(1,:)),'-k')
plot(t,imag(M_t(1,:)),'-r')