function [y]=fun(v,x)
tau=v(1);
tau2=v(2);
y=exp(-x./tau).*exp(-x./tau2);