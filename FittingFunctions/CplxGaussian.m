function [y]=CplxGaussian(k,x)
phi=k(1);
sig=k(2);
A=k(3);
t0=k(4);
y=A*exp(-2*pi*1i*phi.*(x-t0)).*exp(-(x-t0).^2./(2*sig^2));