function [xfit,ypred,beta,pm,resid] = monodecay_t2fit(xdata,ydata,guesses,CI,opts)

% fitting routine and plotting
[beta,resid,J] = nlinfit(xdata,ydata,@t2monofit,guesses,opts);

%output: beta is coefficients
% resid is residuals
% J is Jacobian, COVB is covariance matrix, mse is mean squared error

alpha = 1 - CI/100;

ci = nlparci(beta,resid,'jacobian',J,'alpha',alpha);

% [ypred,delta] = nlpredci(@T1_recovery,x_fit,beta,resid,J);

xfit = 0:max(xdata)/1000:max(xdata);
ypred = nlpredci(@t2monofit,xfit,beta,resid,J);
pm = abs(ci(:,2)-beta(:,1));
% pm(:,1) = abs(ci(:,2) - beta(:,1));

% beta = [beta, pm];




