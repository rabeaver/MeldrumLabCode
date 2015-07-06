function [y]=fun(v,x)
%y0=v(1);
A=v(1);
tau=v(2);
y=(A).*(exp(-x./tau));