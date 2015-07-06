function [xfit,ypred,beta,pm,resid] = monodecay_t2fit_simple(xdata,ydata,guesses,CI); %,opts)

% fitting routine and plotting
[beta,resid,J] = nlinfit(xdata,ydata,@simpleexpdecay,guesses); %,opts);

%output: beta is coefficients
% resid is residuals
% J is Jacobian, COVB is covariance matrix, mse is mean squared error

alpha = 1 - CI/100;

ci = nlparci(beta,resid,'jacobian',J,'alpha',alpha);

% [ypred,delta] = nlpredci(@T1_recovery,x_fit,beta,resid,J);

xfit = 0:max(xdata)/1000:max(xdata);
ypred = nlpredci(@simpleexpdecay,xfit,beta,resid,J);

pm(:,1) = abs(ci(:,2) - beta(:,1)); 





