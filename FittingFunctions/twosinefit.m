function [y]=twosinefit(v,x)
ph1=v(1);
freq1=v(2);
A1=v(3);
tau1=v(4);

ph2=v(5);
freq2=v(6);
A2=v(7);
tau2=v(8);
y=exp(-x./tau1).*A1.*sin(freq1.*x+ph1) + exp(-x./tau2).*A2.*sin(freq2.*x+ph2);