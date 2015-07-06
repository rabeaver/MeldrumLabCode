function [xfit,ypred,beta,resid] = bidecay_t2fit_bound(xdata,ydata,guesses,bounds,CI)

% function [xfit,ypred,beta,resid] = bidecay_t2fit_bound(xdata,ydata,guesses,CI)
% fitting routine and plotting
[x,~,resid,~,~,~,J] = lsqcurvefit(@t2bifit,guesses,xdata,ydata,bounds(:,1),bounds(:,2));

%output: beta is coefficients
% resid is residuals
% J is Jacobian, COVB is covariance matrix, mse is mean squared error

alpha = 1 - CI/100;

ci = nlparci(x,resid,'jacobian',J,'alpha',alpha);

% [ypred,delta] = nlpredci(@T1_recovery,x_fit,beta,resid,J);

xfit = 0:max(xdata)/1000:max(xdata);
ypred = nlpredci(@t2bifit,xfit,x,resid,J);

pm(:,1) = abs(ci(:,2) - x(:,1)); 

beta = [x, pm];




