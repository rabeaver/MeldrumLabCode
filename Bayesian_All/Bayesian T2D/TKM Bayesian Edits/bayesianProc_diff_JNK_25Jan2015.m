% Demonstration of the Bayesian estimation algorithm for multi-component 
% T2 estimation.
%
% Author: Kelvin Layton (klayton@unimelb.edu.au)
% Date: June 2012
% Modified by T Meldrum, July 2014
% 23 Feb 2015, added T2-D plotting
clear
close all
clc

try
    parpool
catch
end

% User parameters
%
datafile= '/Users/tyler/Dropbox/Coding/T2D Solving/Bayesian Recent Edits/Simulated Data/dataSimT2D.mat';

%params for simulated data
T2sim = 0.005; %s
Dsim = 3.3e-9; %m^2 s^-1

%guesses for T2-D
T2guess = [0.002]; %s
Dguess = [3.1e-09];         % diffusion guess (m^2/s)

%params for both (known, measured, no uncertainty)
DELTA = 2e-3;   %s
deltaMin = 10e-6; %s
deltaMax = 710e-6; %s
nEchoes = 32;
nDPnts = 9;
tEcho   = 150e-6;   % echo time (s)
SNR     = 50;      % approximate SNR
nTrials = 50;        % times to iterate noise
T2rangeMin = 1e-4;    % minimum T2 value (s)
T2rangeMax = 1e-1;    % maximum T2 value (s)
DrangeMin = 1e-10;    % minimum D value (m2/s) 
DrangeMax = 1e-8;     % minimum D value (m2/s)
gridCt = 1000;           % number of points in each direction of T2/D
G = 7;         % field gradient T/m
gamma = 267.513e6; % Gyromagnetic Ratio (rad T-1 S-1)

deltaVec = linspace(deltaMin, deltaMax, nDPnts); %s
decayIn = exp(-gamma.^2 .* G.^2 .* deltaVec.^2 .* Dguess .* (DELTA + (deltaVec./3))); % Initial decay


[spec] = T2Dsim(nEchoes, tEcho, T2sim, G, DELTA, Dsim, deltaMin, deltaMax, nDPnts, SNR);

% Read and set up data file
% pw = 6;
% d4 = 40;
% [aq, spec] = tecmag_t2d_1acqu(datafile, nEchoes, tEcho*10^6, pw, d4);

% load(datafile);

specRe = real(spec);
specIm = randn(size(spec))/SNR;
spec = reshape([specRe, specIm]', 2*nDPnts*nEchoes, 1);

s = spec./max(spec);
% plot(real(s));

%%

% Define prior (in terms of R1 and R2)
if length(T2guess) == 1
    %[weight, 1/T2guess, Dweight, Dguess, alpha, phase] 
    opts.priorMean = [1; 1/T2guess; Dguess; pi/2; 0];
elseif length(T2guess) == 2
    %[weight_1, weight_2, 1/T2guess_1, 1/T2guess_2 alpha, phase]
    opts.priorMean = [0.5; 0.5; 1/T2guess(1); 1/T2guess(2); Dguess(1); Dguess(2); pi/2; 0];
end 
opts.priorCov = diag([1^2*ones(length(T2guess),1); 10e4^2*ones(length(T2guess),1); 10e-2^2*ones(length(Dguess),1); 100^2*ones(2,1)]);
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
gammas = logspace(-8,1,150);
gammas = gammas ./ sum(gammas);
opts.correctionSchedule = gammas;

% Allocate output arrays
T2Hat = zeros(length(T2guess),nTrials);
wHat = zeros(length(T2guess),nTrials);
DHat = zeros(length(Dguess), nTrials);
iTrial = 1;

%%
% Repeat for different noise realisations

parfor iTrial=1:nTrials
    
    % Add noise to signal
    y = s + (1./SNR)*randn(size(s));

    % Estimate distribution components
    [T2Hat(:,iTrial), wHat(:,iTrial), DHat(:,iTrial)] = bayesian_estimate(y,opts);
    
    % Keep H (from epg) for each different set of added noise
%     Hout(:,iTrial) = reshape(H,nEchoes*length(T2guess),1);

end


% if length(T2guess) == 2
%     Hout = reshape(Hout, nEchoes,length(T2guess),nTrials);
% end



%%
% Calculate normal distributions
T2grid = logspace(log10(T2rangeMin),log10(T2rangeMax),gridCt);
DGrid  = logspace(log10(DrangeMin),log10(DrangeMax),gridCt);
T2mean = nanmean(T2Hat');
DMean  = nanmean(DHat');
T2std  = nanstd(T2Hat');
DStd   = nanstd(DHat');
if length(T2guess) == 1
    G = DMean*exp(-((T2grid-T2mean).^2)/(T2std/(2*sqrt(2*log(2))))^2);
    GD = T2mean*exp(-((DGrid-DMean).^2)/(DStd/(2*sqrt(2*log(2))))^2);
elseif length(T2guess) == 2
    G = DMean(1)*exp(-((T2grid-T2mean(1)).^2)/(T2std(1)/(2*sqrt(2*log(2))))^2) + ...
        DMean(2)*exp(-((T2grid-T2mean(2)).^2)/(T2std(2)/(2*sqrt(2*log(2))))^2);
    GD = T2mean(1)*exp(-((DGrid-DMean(1)).^2)/(DStd(1)/(2*sqrt(2*log(2))))^2) + ...
         T2mean(2)*exp(-((DGrid-DMean(2)).^2)/(DStd(2)/(2*sqrt(2*log(2))))^2);
end

% Plot estimates
figure(1)
hold on
if length(T2guess) == 1
    plot(T2Hat(1,:),DHat(1,:),'bx');
elseif length(T2guess) == 2
    plot(T2Hat(1,:),DHat(1,:),'bx');
    plot(T2Hat(2,:),DHat(2,:),'rx'); 
end
% plot(T2grid,G,'-k','LineWidth',2);
% plot(Gw,DGrid,'-k','LineWidth',2);
% stem(T2mean,DMean,'Color',0.5*[1 1 1],'LineWidth',2);

xlabel('T_2 [s]')
xlim([T2rangeMin T2rangeMax])
ylim([DrangeMin DrangeMax])
set(gca,'yscale','log','xscale','log')
ylabel('D [m^2/s]')
title('T_2-D')

%
% t = (1:nEchoes)*tEcho;
% plot(t,s)
% cftool(t,s)