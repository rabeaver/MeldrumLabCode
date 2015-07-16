function [y]=fun(k,x)
A=k(1);
t2=k(2);
y=A.*(exp(-x./t2));