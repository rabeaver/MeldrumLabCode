function [y]=fun(k,x)
y0=k(1);
A=k(2);
t2=k(3);
y=y0 + (A).*(exp(-x./t2));
 %y=(A).*(exp(-x./t2));