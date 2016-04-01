function [b,u] = house(column,pivot,systemmat,tailledata,u)
% rotates a vector from "mat" (indexed by column)
% the outputs "u" (vector) and "b" (float) allow us to apply orthogonal tranform elsewhere

% som = 0;
u(1:pivot-1)=0;
u(pivot+1:tailledata)=systemmat(pivot+1:tailledata,column);
us = u.*u;
som = sum(us(pivot+1:tailledata));
if systemmat(pivot,column)<0
    s = sqrt(som + (systemmat(pivot,column))^2);
else
    s = - sqrt(som + (systemmat(pivot,column))^2);
end
u(pivot) = systemmat(pivot,column)-s;
b = s*u(pivot);