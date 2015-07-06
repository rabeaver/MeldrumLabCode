function [y]=fun(v,x)
c=v(1);
% A = v(2);
y=exp(-c*x.^2);