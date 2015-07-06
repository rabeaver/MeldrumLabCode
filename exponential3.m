function [y]=fun(v,x)
tau=v(1);
y=(exp(-x./tau));