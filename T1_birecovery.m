function [y]=fun(k,x)
A1=k(1);
t1_1=k(2);
A2 = k(3);
t1_2 = k(4);
y = A1*(1-exp(-x./t1_1)) + A2*(1-exp(-x./t1_2));