function [y]=fun(k,x)
y0=k(1);
A=k(2);
t1=k(3);
y=y0 + A*(1-exp(-x./t1));