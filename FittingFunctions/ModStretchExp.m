function [y]=fun(k,x)
A    = k(1);
t2   = k(2);
tau  = k(3);
beta = k(4);
y=A*exp(-(x./t2).*((1+(x./tau)).^(beta-1)));