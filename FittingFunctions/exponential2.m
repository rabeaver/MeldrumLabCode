function [y]=fun(v,x)
y0=v(1);
tau=v(2);
y=y0 + (exp(-x./tau));