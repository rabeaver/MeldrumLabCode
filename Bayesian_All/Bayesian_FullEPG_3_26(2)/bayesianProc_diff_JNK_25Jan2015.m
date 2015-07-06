% Demonstration of the Bayesian estimation algorithm for multi-component 
% T2 estimation.
%
% Author: Kelvin Layton (klayton@unimelb.edu.au)
% Date: June 2012
% Modified by T Meldrum, July 2014
%
clear
close all
clc


% Need to figure out how to use the covariance to my advantage. So far, it
% seems to get really crazy unless the covariance scalar factor before the
% D term's covariance is set low (e.g. 1e-4). When it's that low, Psi(3,3)
% doesn't deviate from that value, which means that the D constant never
% deviates and, consequently, the T2 doesn't deviate. It may be just a
% matter of getting a better number, but I fear that it's related to the
% relative effects of diffusion (rapid signal attenuation during a fixed
% time interval) and T2 (relative little signal attenuation during a fixed
% time interval). This *might* be fixed by scaling the gradient, but that
% would reduce the total amount of diffusion, not just increase the
% diffusion constant D. Not sure how to make this work...

% 24 March 2015 The problem of badly scaled matrices arises when the
% estimated value for D (from mu or d_prev) becomes negative. THe next
% matrix calculation gives the badly scaled error. How to avoid this?


% 
% try
%     parpool
% catch
% end

% User parameters
%
datafile= '/Users/tyler/Dropbox/Coding/T2D Solving/Bayesian Recent Edits/Simulated Data/dataSimT2D.mat';

%params for simulated data
T2sim = 0.005; %s
Dsim = 3.3e-9; %m^2 s^-1

%guesses for T2-D
T2guess = [0.005]; %s
Dguess = [3.3e-9];         % diffusion guess (m^2/s)

%params for both (known, measured, no uncertainty)
DELTA = 2e-3;   %s
deltaMin = 15e-6; %s
deltaMax = 815e-6; %s
nEchoes = 32;
nDPnts = 11;
tEcho   = 150e-6;   % echo time (s)
SNR     = 20;      % approximate SNR
nTrials = 10;        % times to iterate noise
T2rangeMin = 10e-4;    % minimum examination range (s)
T2rangeMax = 10e-3;    % maximum examination range (s)
DrangeMin  = 1e-09;
DrangeMax  = 1e-08;
G = 7;         % field gradient T/km
gamma = 267.513e6; % Gyromagnetic Ratio (rad T-1 S-1)

deltaVec = logspace(log10(deltaMin), log10(deltaMax), nDPnts); %s
decayIn = exp(-gamma.^2 .* G.^2 .* deltaVec.^2 .* Dguess .* (DELTA + (2*deltaVec./3))); % Initial decay


[spec] = T2Dsim(nEchoes, tEcho, T2sim, G, DELTA, Dsim, deltaMin, deltaMax, nDPnts, SNR);
[specguess] = T2Dsim(nEchoes, tEcho, T2guess, G, DELTA, Dguess, deltaMin, deltaMax, nDPnts, SNR);

% figure(3)
% hold on
% plot(spec(:,1))
% plot(specguess(:,1))
% legend('spec','specguess')

%
% Read and set up data file
% pw = 6;
% d4 = 40;
% [aq, spec] = tecmag_t2d_1acqu(datafile, nEchoes, tEcho*10^6, pw, d4);

% load(datafile);

specRe = real(spec);
specIm = randn(size(spec))/SNR;
spec = reshape([specRe, specIm]', 2*nDPnts*nEchoes, 1);

specReGuess = real(specguess);
specImGuess = specIm;
specguess = reshape([specReGuess, specImGuess]', 2*nDPnts*nEchoes, 1);

s = spec./max(spec);
sguess = specguess./max(specguess);
figure(1)
hold on
plot(real(s));
plot(real(sguess),'-r')
legend('Simuated','Guess')

%%

% Define prior (in terms of R1 and R2)
if length(T2guess) == 1
    %[weight, 1/T2guess, Dweight, Dguess, alpha, phase] 
    opts.priorMean = [1; 1/T2guess; Dguess; pi/2; 0];
elseif length(T2guess) == 2
    %[weight_1, weight_2, 1/T2guess_1, 1/T2guess_2 alpha, phase]
    opts.priorMean = [0.5; 0.5; 1/T2guess(1); 1/T2guess(2); Dguess(1); Dguess(2); pi/2; 0];
end 

opts.priorCov = diag([10^2*ones(length(T2guess),1); 10000^2*ones(length(T2guess),1); 1000^2*ones(length(Dguess),1); 10^2*ones(2,1)]);
opts.tau = tEcho;
opts.sigma2 = (1./SNR)^2;
opts.D = Dguess;
opts.G = G;
opts.nEchoes = nEchoes;
opts.nDPnts = nDPnts;
opts.decayIn = decayIn;
opts.gamma = gamma;
opts.deltaVec = deltaVec;
opts.DELTA = DELTA;


% Define correction schedule
gammas = logspace(-6,5,70);
gammas = gammas ./ sum(gammas);
opts.correctionSchedule = gammas;

% Allocate output arrays
T2Hat = zeros(length(T2guess),nTrials);
wHat = zeros(length(T2guess),nTrials);
DHat = zeros(length(Dguess), nTrials);
iTrial = 1;

%%
% Repeat for different noise realisations

for iTrial=1:nTrials
    
    % Add noise to signal
    y = s + (1./SNR)*randn(size(s));

    % Estimate distribution components
    [T2Hat(:,iTrial), wHat(:,iTrial), DHat(:,iTrial)] = bayesian_estimate(y,opts);
    [spectemp] = T2Dsim(nEchoes, tEcho, T2Hat(:,iTrial), G, DELTA, DHat(:,iTrial), deltaMin, deltaMax, nDPnts, SNR);
    % Keep H (from epg) for each different set of added noise
%     Hout(:,iTrial) = reshape(H,nEchoes*length(T2guess),1);
    specRe = real(spectemp);
    specIm = randn(size(spectemp))/SNR;
    specnew = reshape([specRe, specIm]', 2*nDPnts*nEchoes, 1);
    snew(:,iTrial) = specnew./max(specnew);
    
end


figure(4)
hold on
plot(real(s),'-k','LineWidth',2)
plot(real(sguess),'-r','LineWidth',2)
plot(real(snew))

% if length(T2guess) == 2
%     Hout = reshape(Hout, nEchoes,length(T2guess),nTrials);
% end



%%
% Calculate normal distributions
T2grid = linspace(T2rangeMin,T2rangeMax,10000);
DGrid  = linspace(DrangeMin,DrangeMax,10000);
T2mean = nanmean(T2Hat');
DMean  = nanmean(DHat');
T2std  = nanstd(T2Hat');
DStd   = nanstd(DHat');
if length(T2guess) == 1
    G = DMean*exp(-((T2grid-T2mean).^2)/(T2std/(2*sqrt(2*log(2))))^2);
    Gw = T2mean*exp(-((DGrid-DMean).^2)/(DStd/(2*sqrt(2*log(2))))^2);
elseif length(T2guess) == 2
    G = DMean(1)*exp(-((T2grid-T2mean(1)).^2)/(T2std(1)/(2*sqrt(2*log(2))))^2) + ...
        DMean(2)*exp(-((T2grid-T2mean(2)).^2)/(T2std(2)/(2*sqrt(2*log(2))))^2);
    Gw = T2mean(1)*exp(-((DGrid-DMean(1)).^2)/(DStd(1)/(2*sqrt(2*log(2))))^2) + ...
         T2mean(2)*exp(-((DGrid-DMean(2)).^2)/(DStd(2)/(2*sqrt(2*log(2))))^2);
end

% Plot estimates
figure(2)
hold on
if length(T2guess) == 1
    plot(log10(T2Hat(1,:)),log10(DHat(1,:)),'bx');
    plot(log10(T2sim),log10(Dsim),'ko','MarkerFaceColor','black');
    plot(log10(T2guess),log10(Dguess),'ro');
    plot(log10(T2mean),log10(DMean),'go','MarkerFaceColor','green');
elseif length(T2guess) == 2
    plot(log10(T2Hat(1,:)),log10(DHat(1,:)),'bx');
    plot(log10(T2Hat(2,:)),log10(DHat(2,:)),'rx'); 
end
plot(log10(T2grid),log10(G),'-k','LineWidth',1);
plot(log10(Gw),log10(DGrid),'-k','LineWidth',1);
% stem(T2mean,DMean,'Color',0.5*[1 1 1],'LineWidth',2);

xlabel('log(T2 times [s])')
xlim(log10([T2rangeMin T2rangeMax]))
ylim(log10([DrangeMin DrangeMax]))
ylabel('log(D [m^2 s^{-1}])')
% title('EPG and raw data comparison')


%%
t = (1:nEchoes)*tEcho;
plot(t,s)
cftool(t,s)