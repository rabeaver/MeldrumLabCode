function [y]=fun(v,x)
y0=v(1);
A=v(2);
tau=v(3);
y=y0 + (A).*(exp(-x./tau));