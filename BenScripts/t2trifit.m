function [y]=fun(k,x)
A_1=k(1);
t2_1=k(2);
A_2=k(3);
t2_2=k(4);
A_3=k(5);
t2_3=k(6);
y= (A_1)*(exp(-x./t2_1)) + (A_2).*(exp(-x./t2_2)) + (A_3).*(exp(-x./t2_3));