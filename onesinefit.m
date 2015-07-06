function [y]=onesinefit(v,x)
ph1=v(1);
freq1=v(2);
A1=v(3);
tau1=v(4);
y=exp(-x./tau1).*A1.*sin(freq1.*x+ph1);