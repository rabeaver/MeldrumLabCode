function [y]=fun(k,x)
Sh20=k(1);
A_1=k(2);
t2_1=k(3);
A_2=k(4);
t2_2=k(5);
% A_3=k(6);
% t2_3=k(7);
y= (Sh20)*((0.9587)*(exp(-x./0.0064)) + (0.0904).*(exp(-x./0.0005)))+ (A_1)*(exp(-x./t2_1)) + (A_2).*(exp(-x./t2_2));% + (A_3).*(exp(-x./t2_3));