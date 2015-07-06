% Programm MULTIEXP to produce Multi exponential data 
% or process them by Laplace inversion with a regularization function
clear all;
tempsrelax = exp((-3 : 0.2 : 1)*log(10));
sigma1 = 0.5;
mu1 = 1;
amp1 = 10;
sigma2 = 0.3;
mu2 = 0.02;
amp2 = 3;
alpha = 10^5;
distribution= amp1*(1/(sigma1*sqrt(2*pi))*exp(-(log(tempsrelax)-log(mu1)).^2/(2*sigma1^2)))+amp2*(1/(sigma2*sqrt(2*pi))*exp(-(log(tempsrelax)-log(mu2)).^2/(2*sigma2^2)));
temps = 0.001 : 0.08 : 10;
taille = size(temps);
for m = 1 : taille(2)
signal(m) = sum(distribution.*exp(-tempsrelax.^(-1)*temps(m)));
end

figure
% plot(temps,signal)
ampmax = max(signal);
bruit = ampmax/100*randn(1,taille(2));
signal = signal + bruit;
% hold on
plot(temps,signal,'r')
% hold off
rel = size(tempsrelax);
tp = size(temps);
figure

semilogx(tempsrelax,distribution)

signalplus = zeros(1,rel(2));
signal=[signal signalplus]
Matriceelt = zeros(tp(2)+rel(2),rel(2));

for k = 1 : tp(2)
    for l = 1 : rel(2)
        Matriceelt(k,l) = exp(-temps(k)/tempsrelax(l));
    end;
end;

for k = 1 : rel(2)
    Matriceelt(tp(2)+k,k) = alpha;
end
    
taille=size(Matriceelt);
x0 = zeros(1,taille(2));
% for t =1:10
    x = lsqnonneg(Matriceelt,signal',x0);
    x0=x;
   
figure
semilogx(tempsrelax,distribution)
hold on
    semilogx(tempsrelax,x)
    axis([0.001,1000,min(x),max(x)])
    hold off
% % end;



hold off
