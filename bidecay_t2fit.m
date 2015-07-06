function [xfit,ypred,beta,beta_err,resid] = bidecay_t2fit(xdata,ydata,guesses,CI) %,opts)

% fitting routine and plotting
[beta,resid,J,~,mse] = nlinfit(xdata,ydata,@t2bifit,guesses); %,opts);

%output: beta is coefficients
% resid is residuals
% J is Jacobian, COVB is covariance matrix, mse is mean squared error

alpha = 1 - CI/100;

ci = nlparci(beta,resid,'jacobian',J,'alpha',alpha);

% [ypred,delta] = nlpredci(@T1_recovery,x_fit,beta,resid,J);

xfit = 0:max(xdata)/1000:max(xdata);
ypred = nlpredci(@t2bifit,xfit,beta,resid,J);

pm(:,1) = abs(ci(:,2) - beta'); 
beta_err = pm;
% beta = [beta, pm];




