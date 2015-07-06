function [A,tau,ycalc,resid,xpred,ypred] = linearizedT1fit(xdata,ydata)
warning off all
                                %diasble warnings so the errors of the nlin fit for determining the yoffset won't go crazy

guesses = [ max(ydata);      %initial guesses for determinig the yoffset
            max(ydata);
            xdata(round(length(xdata)/5)) ];
        
xpred = 0:max(xdata)/1000:max(xdata); %a time scale with 1001 points over the same range as the original data for assessing the fit

try
    beta = nlinfit(xdata,ydata,@T1_recovery,guesses);     %try the nlinfit to get a yoffset value
catch
    A = 0;                                              % if it fails, just zero out the data set--it's probably useless
    tau = 0;
    ycalc = zeros(length(xdata),1);
    resid = zeros(length(xdata),1);
    ypred = zeros(length(xpred),1);
%     fprintf('Failed y0 exponential fit.\n')
    return
end
warning on all                  %resume warnings

y0 = beta(1);
logx = log(xdata);           %make a linearized version of the ydata for linear fitting
newY = ydata - y0;


p = polyfit(logx,newY,1);       %fit the data to a line

A = p(2);      %the scaling factor
tau = 1/p(1);      %the decay constant

ycalc = y0 + A*(1-exp(-xdata./tau));   %calculate the data points at the measured times using the fit and the residuals

resid = ydata - ycalc;

figure
hold on
scatter(xdata,ydata)
scatter(xdata,ycalc,'or')

ypred = y0 + A*(1-exp(-xpred./tau));   %calculate the data using the 1001 points
end