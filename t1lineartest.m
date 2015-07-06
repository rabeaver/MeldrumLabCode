close all
for n = 1:1:8

figure(n)
subplot(2,1,1)
scatter(T1time,real(allT1Phased(:,n)))
subplot(2,1,2)
scatter(log(T1time),real(allT1Phased(:,n)))

p = polyfit(log(T1time),real(allT1Phased(:,n)),1);

tau(n) = 1/p(1);
A(n) = p(2);

xfit = 0:max(T1time)/1000:max(T1time);
yfit = A(n)*(1-exp(-xfit./tau(n)));

figure(n)
subplot(2,1,1)
hold on
% scatter(T1time,real(allT1Phased(:,n)))
plot(xfit,yfit)
end