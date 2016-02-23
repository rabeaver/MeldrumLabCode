function [xfit,yfit,coeffs,resid,pm] = T1_fit_no_offset(xdata,ydata,guesses,CI)

% fitting routine and plotting
[beta,resid,J,COVB,mse] = nlinfit(xdata,ydata,@T1_recovery_no_offset,guesses);

%output: beta is coefficients
% resid is residuals
% J is Jacobian, COVB is covariance matrix, mse is mean squared error

alpha = 1 - CI/100;
ci = nlparci(beta,resid,'jacobian',J,'alpha',alpha);

% [ypred,delta] = nlpredci(@T1_recovery,x_fit,beta,resid,J);

xfit = 0:max(xdata)/1000:max(xdata);
yfit = beta(1) - beta(1).*exp(-xfit./beta(2));

% S10/80
% [ypred1080,delta1080] = nlpredci(@T1_recovery_no_offset,[10;80],beta,resid,J,alpha);
% s1080(1) = ypred1080(1)/ypred1080(2);
% s1080(2) = sqrt((delta1080(1)/ypred1080(2))^2 - (ypred1080(1)*delta1080(2)/(ypred1080(2)^2))^2);

coeffs(:,1) = beta;
pm = ci(:,2)-beta;
coeffs(:,2) = pm;



