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

try
    parpool
catch
end

% User parameters
%
datafile= '/Users/jaredking/Documents/Research Files and Data/NewBayesianStuff_12_06_2014/Simulated Data/dataSimT2D.mat';


T2guess = [0.02]; %s
Dguess = [3e-09];         % diffusion guess (m^2/s)
DELTA = 2e-3;   %s
deltaVec = linspace(10e-6, 310e-6, 21); %s
nEchoes = 64;
nDPnts = 21;
tEcho   = 120e-6;   % echo time (s)
SNR     = 500;      % approximate SNR
nTrials = 7;        % times to iterate noise
rangeMin = 10e-6;    % minimum examination range (s)
rangeMax = 110e-6;    % maximum examination range (s)
G = 7;         % field gradient T/m
gamma = 267.513e6; % Gyromagnetic Ratio (rad T-1 S-1)
decayIn = exp(-gamma.^2 .* G.^2 .* deltaVec.^2 .* Dguess .* (DELTA + (deltaVec./3))); % Initial decay

% Read and set up data file
pw = 6;
d4 = 40;
[aq, spec] = tecmag_t2d_1acqu(datafile, nEchoes, tEcho*10^6, pw, d4);
specRe = real(spec);
specIm = randn(size(spec))/20;
spec = reshape([specRe, specIm]', 2*nDPnts*nEchoes, 1);

s = spec./max(spec);
plot(real(s));

% Define prior (in terms of R1 and R2)
if length(T2guess) == 1
    %[weight, 1/T2guess, Dweight, Dguess, alpha, phase] 
    opts.priorMean = [1; 1/T2guess; Dguess; pi/2; 0];
elseif length(T2guess) == 2
    %[weight_1, weight_2, 1/T2guess_1, 1/T2guess_2 alpha, phase]
    opts.priorMean = [0.5; 0.5; 1/T2guess(1); 1/T2guess(2); Dguess(1); Dguess(2); pi/2; 0];
end 
opts.priorCov = diag([1^2*ones(length(T2guess),1); 10000^2*ones(length(T2guess),1); (1e-03)^2*ones(length(Dguess),1); 100^2*ones(2,1)]);
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
gammas = logspace(-6,1,30);
gammas = gammas ./ sum(gammas);
opts.correctionSchedule = gammas;

% Allocate output arrays
T2Hat = zeros(length(T2guess),nTrials);
wHat = zeros(length(T2guess),nTrials);
DHat = zeros(length(Dguess), nTrials);
iTrial = 1;

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
T2grid = linspace(rangeMin,rangeMax,1000);
wGrid  = linspace(0,1.5,1000);
T2mean = mean(T2Hat');
wMean  = mean(wHat');
T2std  = std(T2Hat');
wStd   = std(wHat');
if length(T2guess) == 1
    G = wMean*exp(-((T2grid-T2mean).^2)/(T2std/(2*sqrt(2*log(2))))^2);
    Gw = T2mean*exp(-((wGrid-wMean).^2)/(wStd/(2*sqrt(2*log(2))))^2);
elseif length(T2guess) == 2
    G = wMean(1)*exp(-((T2grid-T2mean(1)).^2)/(T2std(1)/(2*sqrt(2*log(2))))^2) + ...
        wMean(2)*exp(-((T2grid-T2mean(2)).^2)/(T2std(2)/(2*sqrt(2*log(2))))^2);
    Gw = T2mean(1)*exp(-((wGrid-wMean(1)).^2)/(wStd(1)/(2*sqrt(2*log(2))))^2) + ...
         T2mean(2)*exp(-((wGrid-wMean(2)).^2)/(wStd(2)/(2*sqrt(2*log(2))))^2);
end

% Plot estimates
figure(1)
hold on
if length(T2guess) == 1
    plot(T2Hat(1,:),wHat(1,:),'bx');
elseif length(T2guess) == 2
    plot(T2Hat(1,:),wHat(1,:),'bx');
    plot(T2Hat(2,:),wHat(2,:),'rx'); 
end
plot(T2grid,G,'-k','LineWidth',2);
plot(Gw,wGrid,'-k','LineWidth',2);
stem(T2mean,wMean,'Color',0.5*[1 1 1],'LineWidth',2);

xlabel('T2 times [s]')
xlim([rangeMin rangeMax])
ylabel('weight [arb]')
title('EPG and raw data comparison')


%%
t = (1:nEchoes)*tEcho;
plot(t,s)
cftool(t,s)