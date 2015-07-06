close all
clc
clear

%% Load data
if ispc == 1;
    filedir = 'C:\Users\jnking01\Desktop\'; 
    filename = 'C:\Users\jnking01\Desktop\MilkT2D_2.tnt';
elseif ismac == 1;
    filedir = '/Users/tyler/Dropbox/Data/T2D_Water_T2D/';
    filename = '/Users/tyler/Dropbox/Data/T2D_Water_T2D/H2O_T2D_test3_10Oct2014.tnt';
end

cd(filedir)

[aq,spec,spec2] = readTecmag4d(filename);

n1D = aq.td(1);
n2D = aq.td(2);
nPts = aq.ctd;
nEchoes = 1024;
tEcho = 150e-6;     % echo time
tD = 1e-6;
ptsPer_tE = (tEcho/tD);
ptsPerEcho = 69;

dMin = 45e-6; %minimum delta time in s
DELTA = 1e-3; %DELTA in s
dStep = 27e-6; %step size between deltas in s

G = 6.5998;         % T/m (Smallest gradient amplitude) (6.5998 PM25; 23.8626 PM5)


tauVecA = dMin:dStep:dMin+(dStep*(n2D-1));
tauVec = tauVecA;

spec2 = [spec2,zeros(n2D,int32(nEchoes*ptsPer_tE-size(spec2,2)))];
spec3 = reshape(spec2',int32(ptsPer_tE),nEchoes,n2D);
spec3 = spec3(1:ptsPerEcho,:,:);

dataInt = reshape(sum(spec3,1),nEchoes*n2D,1);
dataIntSurf = reshape(sum(spec3,1),nEchoes,n2D);

dISRe = real(dataIntSurf);
dISIm = imag(dataIntSurf);

dISRe2 = reshape(dISRe,size(dISRe,1)*size(dISRe,2),1);
dISIm2 = reshape(dISIm,size(dISIm,1)*size(dISIm,2),1);

dIS = [dISRe2';dISIm2'];

dISnew = reshape(dIS,2*nEchoes,n2D);
dISnew = dISnew';


save('dataIntSurf7.dat','dISnew','-ascii','-tabs');
fitData = real(dataInt)./max(real(dataInt));
% fitData = real(dataIntSurf(:,1))./max(real(dataIntSurf(:,1)));

echoVector = (1:nEchoes)*tEcho;

figure(2)
hold on
plot(fitData,'-k')

%% Calc data for 1st set
% if ismac == 1;
%     epgdir = '/Users/tyler/Dropbox/Coding/Nick_Coding/epg_funcs/';
% elseif ispc ==1;
%     epgdir = 'C:\Users\tkmeldrum\Dropbox\Coding\Nick_Coding\epg_funcs\';
% end

% cd(epgdir)

D = 2.2e-9;         % Diffusion (m^2/s)
T1 = inf;       % s
T2 = 100e-3;         % s
ang1 = pi/2;
% n2D = 1;

k = [D, T2];
% k = [T1];
p = [DELTA,tEcho,n2D,nEchoes,G,T1,ang1];

% ii = 3;
% data1 = real(dataIntSurf(:,ii))./max(real(dataIntSurf(:,ii)));

% figure(1)
% plot(dataInt)


[S,dK]=EPGDiffCalc(k,p,tauVec);
% dk = gamma*G*tEcho; % dk (gradient winding per tEcho
% dk1 = gamma.*G.*tauVec; %gradient winding during delta
% dk2 = gamma*G*DELTA; %gradient winding during DELTA
% clear('S');

% ====== START SEQUENCE =========
% for j = 1:n2D %loop for tau points
% %     for pp = 1:length(pp1); %loop for phase cycle
%         F = [0;0;1];	% Equilibrium Magnetization.
%         % [F+; F-; Z],  all longitudinal in Z0 state.
%         F = epg_rf(F,pi/2,0);			% 90 RF Rotation.
%         [F,EE,BV] = epg_grelax(F,T1,T2,tauVec(j),dk1(j),D,1,0);	% Relaxation, Diffusion, Gradient
%         F = epg_rf(F,pi/2,0);			% 90 RF Rotation.
%         [F,EE,BV] = epg_grelax(F,T1,T2,DELTA,dk2,D,1,0);	% Relaxation, Diffusion, Gradient
%         F = epg_rf(F,pi/2,0);			% 90 RF Rotation.
%         [F,EE,BV] = epg_grelax(F,T1,T2,tauVec(j),dk1(j),D,1,0);	% Relaxation, Diffusion, Gradient
%         S(1,j) = F(1,1);			% Signal is F+(1) state.
%         
%         for n = 2:nEchoes	% Loop with multiple echoes.
%             [F,EE,BV] = epg_grelax(F,T1,T2,tEcho/2,dk/2,D,1,0);	% Relaxation, Diffusion, Gradient
%             F = epg_rf(F,pi,3*pi/2);			% 180 RF Rotation.
%             [F,EE,BV] = epg_grelax(F,T1,T2,tEcho/2,dk/2,D,1,0);	% Relaxation, Diffusion, Gradient
%             S(n,j) = F(1,1);			% Signal is F+(1) state.
%             
%             % -- Calculate B value
%             % b_val(n) = ((gamma*n*G*tEcho)^2)*(2*tEcho/3);
%         end;
% % end
% end



% S = reshape(sum(S,2),nEchoes,n2D);
% S = reshape(S,nEchoes*n2D,1);
% S = S./max(S);

ls = EPGls(k,tauVec,p,fitData);

figure(2)
hold on
plot(real(S))
plot(ls)



%%

opts = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt'); %fitting options

opts = optimoptions(@lsqnonlin,opts); %fitting options

beta = lsqnonlin('EPGls',k,[], [], opts, tauVec, p, fitData);

newFit = EPGDiffCalc(beta,p,tauVec);

figure(3)
hold on
plot(real(newFit))
plot(fitData)

% phat(:,2) = abs(phat(:,2));
% meanphat = mean(phat,1); %avg value for p-hat (all echo shape parameters)
% meanpsi = meanphat(1);
% stdpsi = std(phat(:,1));
% meansigma = meanphat(2);
% stdsigma = std(phat(:,2));
%%
% figure(3)
% surf(abs(S(:,:,1)))
% shading flat
% 
% figure(4)
% surf(abs(dataIntSurf))
% shading flat

%%
Slin = reshape(S,nEchoes*n2D,1);
Slin = Slin./max(abs(Slin));

figure(5)
hold on
plot(abs(Slin),'-k')
plot(abs(dataInt)./max(abs(dataInt)),'-b')
