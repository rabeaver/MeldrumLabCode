clear
clc
close all

N = 256;
tau = 150e-6;
R1 = 0;
R2vec = 1/0.03;
alpha = pi/2;
D = 2.2e-9;         % Diffusion (m^2/s)
G = 6.5998;         % T/m (Smallest gradient amplitude) (6.5998 PM25; 23.8626 PM5)

t = (1:N)*tau;

[ H0 ] = epg( N, tau, R1, R2vec, alpha, D, G );
[ Jr0, Ja0 ] = epg_derivatives( N, tau, R1, R2vec, alpha , D, G);

D = 1.8e-9;         % Diffusion (m^2/s)

[ HD ] = epg( N, tau, R1, R2vec, alpha, D, G );
[ JrD, JaD ] = epg_derivatives( N, tau, R1, R2vec, alpha , D, G);

figure(1)
hold on
plot(t,H0,'-k')
plot(t,HD,'-r')
legend('without diff','with diff')

figure(2)
subplot(1,2,1)
hold on
plot(t,Jr0,'--k')
plot(t,JrD,'--r')
legend('Jr without','Jr with')
subplot(1,2,2)
hold on
plot(t,Ja0,'-k')
plot(t,JaD,'-r')
legend('Ja with','Ja without')
