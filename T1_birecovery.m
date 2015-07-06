function [y]=fun(k,x)
y0=k(1);
A_1=k(2);
t1_1=k(3);
A_2=k(4);
t1_2=k(5);
y=y0 - (A_1).*(exp(-x./t1_1)) - (A_2).*(exp(-x./t1_2));