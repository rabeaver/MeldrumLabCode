function [y]=fun(k,x)
A=k(1);
t2=k(2);
x0=k(3);
y=(A).*(exp(-(x-x0)./t2));
% y=(A).*(exp(-x./t2));