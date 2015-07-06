function [y]=fun(v,x)
A=v(1);
y0=v(2);
b=v(3);
c=v(4);
y=y0 + A.*(1./(1+exp(-b.*(x-c))));

% Fits to a sigmoidal curve of the form f(x) = y0 + A(1+exp(-b(x-c)))^(-1)