function [xfit,ypred,beta,resid] = tridecay_t2fit(xdata,ydata,guesses,CI)

% fitting routine and plotting
[beta,resid,J] = nlinfit(xdata,ydata,@t2trifit,guesses);

%output: beta is coefficients
% resid is residuals
% J is Jacobian, COVB is covariance matrix, mse is mean squared error

alpha = 1 - CI/100;

ci = nlparci(beta,resid,'jacobian',J,'alpha',alpha);

% [ypred,delta] = nlpredci(@T1_recovery,x_fit,beta,resid,J);

xfit = 0:max(xdata)/1000:max(xdata);
ypred = nlpredci(@t2trifit,xfit,beta,resid,J);

pm(:,1) = abs(ci(:,2) - beta(:,1)); 

beta = [beta, pm];
nrElem = size(beta,1) * size(beta,2);

beta = reshape(beta.',nrElem,1);




