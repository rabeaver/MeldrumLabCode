function [y]=fun(k,x)
A_1=k(1);
t2_1=0.0387;
A_2=k(2);
t2_2=0.0057;
y=(A_1)*(exp(-x./t2_1)) + (A_2).*(exp(-x./t2_2));