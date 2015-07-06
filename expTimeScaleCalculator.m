close all
guess=2.5;
taumin=0.012;
nrPnts = 41;
taufrac = .52;

n=log(guess-taumin+1);
l=taufrac*(nrPnts-1);


k=0:nrPnts-1;
tauaxis = exp(n.*(1/l).*k)-1+taumin;

figure(1)
plot(tauaxis,0,'ok')