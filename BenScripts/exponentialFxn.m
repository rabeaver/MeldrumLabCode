function [y] = exponentialFxn(k,x)
a = k(1);
b = k(2);

y = a*exp(b*x);
end

