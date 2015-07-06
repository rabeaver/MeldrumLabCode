clear
clc
close all

dil = 1000;

A = [44331,   2327.6; %extinction coefficients for BSA, Nap [e(280,BSA) e(280,Nap);
     1136.1,  1744];  %                                      e(331,BSA) e(331,Nap)];
 
B(:,1) = [1.9724; 1.0004]; %measured absorbances [A280; A331] 1:9
B(:,2) = [1.4538; 0.7174]; %measured absorbances [A280; A331] 3:7
B(:,3) = [1.3962; 0.6913]; %measured absorbances [A280; A331] 4:6
B(:,4) = [1.0607; 0.5242]; %measured absorbances [A280; A331] 5:5
B(:,5) = [0.86975;0.41792]; %measured absorbances [A280; A331]6:4
B(:,6) = [0.69865;0.3284]; %measured absorbances [A280; A331] 7:3
B(:,7) = [0.34909;0.1357]; %measured absorbances [A280; A331] 9:1

for n = 1:7;

x = A\B(:,n);
c(:,n) = x*1000*dil;   %[[BSA]; [Nap]] in mM, scaled by the dilution factor above
end

xaxis = [0.1 0.3 0.4 0.5 0.6 0.7 0.9];

figure
scatter(xaxis,c(1,:))
title('BSA')
figure
scatter(xaxis,c(2,:))
title('Nap')
figure
scatter(xaxis,c(1,:)./c(2,:))