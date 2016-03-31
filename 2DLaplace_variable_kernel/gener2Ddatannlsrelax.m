clear all;
tic

% generation of the simulated 2D data set
% first: generation of the 2D distribution and then from this distribution: generation of the multiexponentiel 2D data set
sig1a = 0.1;
sig2a = 0.1;
mu1a = 0.5;
mu2a = 0.5;
ampa = 30

sig1b = 0.2;
sig2b = 0.3;
mu1b = 0.2;
mu2b = 0.2;
ampb = 0;

ig1a = 0.9;
sig2a = 1;
mu1a = 0.08;
mu2a = 0.8;
ampa = 30

sig1b = 0.9;
sig2b = 0.8;
mu1b = 0.9;
mu2b = 0.03;
ampb = 40;


T1st = 15;
T2st = 15;
pT1 = 4/(T1st-1);
pT2 = 4/(T2st-1);
T1vl = -3:pT1:1;
T2vl = -3:pT2:1;
T1v = exp(T1vl*log(10));
T2v = exp(T2vl*log(10));
[a,T1s] = size(T1v);
[a,T2s] = size(T2v);

for T1 = 1 : T1s
        distribution(:,T1) = ampa*(1/(sig1a*sqrt(2*pi)))*(1/(sig2a*sqrt(2*pi)))*exp((-(log(T1v(T1))-log(mu1a))^2/(2*(sig1a^2)))-((log(T2v')-log(mu2a)).^2/(2*(sig2a^2)))) + ampb*(1/(sig1b*sqrt(2*pi)))*(1/(sig2b*sqrt(2*pi)))*exp((-(log(T1v(T1))-log(mu1b))^2/(2*(sig1b^2)))-((log(T2v')-log(mu2b)).^2/(2*(sig2b^2))));
end

figure;
if T1s ~= 1
     size(distribution)
     T1st = 32;
     T2st = 32;
     pT1 = 4/(T1st-1);
     pT2 = 4/(T2st-1);
     XI = -3:pT1:1;
     YI = -3:pT2:1;
     distributioni = interp2(T1vl',T2vl,distribution,XI',YI,'cubic');
    surf(T1vl,T2vl,distribution)
    shading interp
    axis([T1vl(1),T1vl(T1st),T2vl(1),T2vl(T2st)])
    set(gcf,'Renderer','zbuffer')
    colorbar
else
    semilogx(T2v,distribution)
end

savefile = ['distribution3' 'm.fig'];
saveas(gcf,savefile)
savefile = ['distribution3' 'm.jpg'];
saveas(gcf,savefile)

figure
surf(T1vl,T2vl,distribution)
shading interp
set(gcf,'Renderer','zbuffer');
axis([T1vl(1),T1vl(T1st),T2vl(1),T2vl(T2st)])
set(gca,'visible','off');
savefile = ['distribution3' 'm.avi'];
mov = avifile(savefile, 'quality',100)
F = getframe(gca);
mov = addframe(mov,F);
for i=1:3:90
    camorbit(-0.3,-2.2,'data',[0 0.5 1])
    drawnow
    F = getframe(gca);
    mov = addframe(mov,F);
end
mov = close(mov);




toc

% data set generation
% nbtm1 = 30;
% nbtm2 = 30;
% pastm1 = 
tm1 = 0.01 : 0.03 : 10;
tm2 = 0.01 : 0.03 : 10;
[a,tm1s] = size(tm1);
[a,tm2s] = size(tm2);
data = zeros(tm2s,tm1s);
for t1 = 1 : tm1s
    for t2 = 1 : tm2s
        data(t2,t1) = sum(sum((distribution*exp(-tm2(t2)./T2v')).*exp(-tm1(t1)./T1v')));
    end
end
data = data';
tm1 = tm1'
tm2 = tm2'
save tm1.test tm1 -ascii
save tm2.test tm2 -ascii
save data.test data -ascii

% save temp.dat A -ascii

ampmax = max(max(data));
bruit = ampmax/200*randn(tm2s,tm1s);
data = data + bruit;
times = 1;

tm1
size(tm1)
tm1i = interp1(1:tm1s,tm1,1:1/times:tm1s); 
[nimportequoi,tm1is] = size(tm1i)

tm2
size(tm2)
tm2i = interp1(1:tm2s,tm2,1:1/times:tm2s); 
[nimportequoi,tm2is] = size(tm2i)

size(data)
data = interp2((1:tm1s)',(1:tm2s),data,(1:1/times:tm1s)',(1:1/times:tm2s),'nearst');
size(data)


% tm2i = tm2;
% tm1i = tm1;
% tm1is = tm1s;
% tm2is = tm2s;

figure;
if T1s ~= 1
    surf(tm1i,tm2i,data)
    axis([tm1i(1),tm1i(tm1is),tm2i(1),tm2i(tm2is)])
    shading interp
    colorbar
else
    plot(tm2,data)
end
toc

% determination of the parameters of the nnlsrelax algorithm
% definition of the min and max relaxation times for the two dimensions "Ta" "Tb"
% that will be defined in the E matrix

Tamin = 0.001;
Tamax = 10;
Tamm = [Tamin Tamax];      % Ta min and max
stepsa = 15;
Tbmin = 0.001;
Tbmax = 10;
Tbmm = [Tbmin Tbmax];
stepsb = 15;
alpha = 1e12;   %10e7
beta = -1;
orient = 'b'; % 'h' for horizontal, 'v' for vertical and 'b' for both directions

if stepsb == 1
    Tbmin = 1;
    Tbmax = 1;
end

timea = tm1i;
timeb = tm2i;
% NNLS with smoothing computation on the 2d data set
[spectrum,taua,taub,chisq,compte]=upnnlsmooth3Dsvd(data,timea,timeb,Tamm,stepsa,Tbmm,stepsb,alpha,beta,orient);
chisq

%elapse time of the computation
toc

taulb = log10(taub);
stb = size(taulb)
taula = log10(taua);
sta = size(taula)
figure;
if stepsb ~= 1
%     T1st = 32;
%     T2st = 32;
%     pT1 = 4/(T1st-1);
%     pT2 = 4/(T2st-1);
%     XI = -3:pT1:1;
%     YI = -3:pT2:1;
%     spectrumi = interp2(taulb',taula,spectrum,XI',YI,'spline');
XI = taula
YI = taulb
spectrumi = spectrum
    surf(XI,YI,spectrumi);
    shading interp
    set(gcf,'Renderer','zbuffer')
    colorbar
    axis([XI(1),XI(sta(2)),YI(1),YI(stb(2))])
else
    semilogx(taua,spectrum,'r')
    hold on
    semilogx(T2v,distribution)
    hold off;
end

