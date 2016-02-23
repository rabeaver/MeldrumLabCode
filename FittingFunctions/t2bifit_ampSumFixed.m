function [y]=fun(k,x,Aa,Ab)
% y0=k(1);
A_1=k(1);
t2_1=k(2);
t2_2=k(3);
% y=y0 + (A_1)*(exp(-x./t2_1)) + (A_2).*(exp(-x./t2_2));
y=(A_1).*(exp(-x./t2_1)) + (1-A_1).*(exp(-x./t2_2));