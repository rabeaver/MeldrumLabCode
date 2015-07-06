function [y]=fun(k,x)
A=k(1);
t1=k(2);
y= + A*(1-exp(-x./t1));