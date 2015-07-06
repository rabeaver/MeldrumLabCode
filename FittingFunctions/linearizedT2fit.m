function [A,tau,ycalc,resid,xpred,ypred] = linearizedT2fit(xdata,ydata)
warning off all%diasble warnings so the errors of the nlin fit for determining the yoffset won't go crazy

guesses = [ 0;      %initial guesses for determinig the yoffset
            max(ydata);
            xdata(round(length(xdata)/7)) ];
        
xpred = 0:max(xdata)/1000:max(xdata); %a time scale with 1001 points over the same range as the original data for assessing the fit

try
    beta = nlinfit(xdata,ydata,@t2monofit,guesses);     %try the nlinfit to get a yoffset value
catch
    A = 0;                                              % if it fails, just zero out the data set--it's probably useless
    tau = 0;
    ycalc = zeros(length(xdata),1);
    resid = zeros(length(xdata),1);
    ypred = zeros(length(xpred),1);
    fprintf('Failed y0 exponential fit.\n')
    return
end
warning on all %resume warnings

ydata = ydata - beta(1);        %shift the ydata by the yoffset

linearized = log(ydata);        %make a linearized version of the ydata for linear fitting

newY = 0;                       %empty the newY data var
skip = 0;                       % a skip variable for shifting the data in case it should be fitted (user input) but the first pt was neg.

if ~isreal(linearized(1))       %if the first point is imaginary (negative before log), plot and ask if it's worth fitting
    figure(8)
    scatter(xdata,ydata)
    choice = menu('Fit this dataset?','Yes','No');
    close(8)
    if choice == 2
        A = 0;                                              % if it fails, just zero out the data set--it's probably useless
        tau = 0;
        ycalc = zeros(length(xdata),1);
        resid = zeros(length(xdata),1);
        ypred = zeros(length(xpred),1);
        return
    
    else
        counter = 1;            %if the user wants the dataset fitted, do the same as in lines 33-41 but with a skip of 1 to dodge the first point
        while isreal(linearized(counter+1))&&counter<length(linearized)
            newY(counter,1) = linearized(counter+1);
            counter = counter+1;
        end
    skip = 1;
    end
end

if ~isreal(linearized)          %this goes through the linearized dataset and stores each value in newY until it encounters the first imaginary linearized value (negative in real units). If there are no such values, it copies over the whole thing.
    counter = 1;
    while isreal(linearized(counter))&&counter<length(linearized)
       newY(counter,1) = linearized(counter);
       counter = counter+1;
    end
else
    newY = linearized;
end

newX = xdata(1+skip:length(newY)+skip); %make the newX data for fitting the same size as the newY data

p = polyfit(newX,newY,1);       %fit the data to a line

A = exp(p(2));      %the scaling factor
tau = -1/p(1);      %the decay constant

calcY = polyval(p,xdata);   %calculate the data points at the measured times using the fit and the residuals
ycalc = exp(calcY);
resid = ydata - ycalc;

ypred = polyval(p,xpred);   %calculate the data using the 1001 points
ypred = exp(ypred);
end