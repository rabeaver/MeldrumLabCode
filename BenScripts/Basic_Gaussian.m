% Variable for Gaussian
gaussianLb = 1;

% Creates Gaussian Matrices??
% g_x = exp(-(t-(max(t)/2)).^2/(2*(1/lb)^2)); % Supresses noise
c = 2*pi*params.pulseLength/1000/(2*sqrt(2*log(2))); %parameters and fit for gaussian filter
g_x = exp(-((t-params.acqTime/2).^2)/(2*c^2));

g_x = repmat(g_x,params.nrEchoes);
g_x = g_x(1:params.nrEchoes,1:params.nrPts);

% *****************************************************************
% So after that just multiply your real and imaginary datas by g_x
% *****************************************************************