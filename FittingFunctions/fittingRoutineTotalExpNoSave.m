function [yData,p] = fittingRoutineTotalExpNoSave(offOrOn,time,confInt)
yOrig=offOrOn; 

y0 = 0.1;
A=max(yOrig);
tau=10;
vStart=[y0,A,tau];
%fprintf('Start: A=%f  tau=%f\n',vStart(1),vStart(2));
yStart=totalExponential(vStart,time);

[vEnd, r, J, COVB, mse] = nlinfit(time,yOrig,@totalExponential,vStart); %fit raw data points to y=A*e^(-x/tau) to find scaling
yEnd=totalExponential(vEnd,time);
%fprintf('End: A=%f  tau=%f\n', vEnd(1),vEnd(2));

ci = nlparci(vEnd, r, 'covar', COVB,'alpha',(1-confInt)); %generates 90% confidence intervals for the fit coefficients

yData = yEnd;

%% final results and plotting
finalFit = [time, yOrig, yEnd];
fitparams = [vEnd(3), abs(vEnd(3) - ci(3,1)); vEnd(1), abs(vEnd(1) - ci(1,1)); vEnd(2), abs(vEnd(2) - ci(2,1))];

% fprintf('End: tau = %f +/- %f\n     y0 = %f +/- %f\n     A = %f +/- %f\n', vEnd(3), abs(vEnd(3) - ci(3,1)), vEnd(1), abs(vEnd(1) - ci(1,1)), vEnd(2), abs(vEnd(2) - ci(2,1)));
fprintf('End: tau = %f +/- %f\n', vEnd(3), abs(vEnd(3) - ci(3,1)));


p = fitparams(1,:);

