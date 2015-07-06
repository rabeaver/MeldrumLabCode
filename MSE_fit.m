function [xfit,ypred,beta,beta_err,resid] = MSE_fit(xdata,ydata,guesses,CI)

% fitting routine and plotting
[beta,resid,J,~,mse] = nlinfit(xdata,ydata,@ModStretchExp,guesses);

%output: beta is coefficients
% resid is residuals
% J is Jacobian, COVB is covariance matrix, mse is mean squared error

alpha = 1 - CI/100;

ci = nlparci(beta,resid,'jacobian',J,'alpha',alpha);

% [ypred,delta] = nlpredci(@T1_recovery,x_fit,beta,resid,J);

xfit = 0:max(xdata)/1000:max(xdata);
ypred = nlpredci(@ModStretchExp,xfit,beta,resid,J);

pm(:,1) = abs(ci(:,2) - beta(:,1)); 
beta_err = pm;
% beta = [beta, pm];




