function [yData,scaling,scalingError,tau,tauError] = FIDfit(dataToFit,time,confInt)
yOrig=dataToFit; 

A=max(yOrig);
tau=10;
vStart=[A,tau];
%fprintf('Start: A=%f  tau=%f\n',vStart(1),vStart(2));
% yStart=exponential1(vStart,time);

[vEnd, r, J, COVB, mse] = nlinfit(time,yOrig,@exponential1,vStart); %fit raw data points to y=A*e^(-x/tau) to find scaling
yEnd=exponential1(vEnd,time);
%fprintf('End: A=%f  tau=%f\n', vEnd(1),vEnd(2));

ci = nlparci(vEnd, r, 'covar', COVB,'alpha',(1-confInt)); %generates 90% confidence intervals for the fit coefficients

scaling = vEnd(1);
scalingError = abs((ci(1,1)-ci(1,2))/2);

tau = vEnd(2);
tauError = abs((ci(2,1)-ci(2,2))/2);

yData = yEnd;


