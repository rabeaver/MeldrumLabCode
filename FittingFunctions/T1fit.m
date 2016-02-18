function [xfit,yfit,coeffs,resid,J,mse,ci,se] = T1fit(xdata,ydata,guesses,CI)

% fitting routine and plotting
[beta,resid,J,COVB,mse] = nlinfit(xdata,ydata,@T1_recovery,guesses);

%output: beta is coefficients
% resid is residuals
% J is Jacobian, COVB is covariance matrix, mse is mean squared error

alpha = 1 - CI/100;
[ci, se] = nlparci(beta,resid,'jacobian',J,'alpha',alpha);

% [ypred,delta] = nlpredci(@T1_recovery,x_fit,beta,resid,J);

xfit = 0:max(xdata)/1000:max(xdata);
yfit = T1_recovery(beta,xfit);

% % S10/80
% [ypred1080,delta1080] = nlpredci(@T1_recovery,[10;80],beta,resid,J,alpha);
% s1080(1) = ypred1080(1)/ypred1080(2);
% s1080(2) = sqrt((delta1080(1)/ypred1080(2))^2 - (ypred1080(1)*delta1080(2)/(ypred1080(2)^2))^2);

coeffs(:,1) = beta;
pm = ci(:,2)-beta;
coeffs(:,2) = pm;



