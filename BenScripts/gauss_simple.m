function [y]=fun(k,x)
A=k(1);
C=k(2);
y = (A).*exp(-(x.^2/(2*C^2)));