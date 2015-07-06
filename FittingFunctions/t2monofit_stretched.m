function [y]=t2monofit(k,x)
A     = k(1);
t2    = k(2);
alpha = k(3);
y=(A).*(exp(-(x./t2).^alpha));