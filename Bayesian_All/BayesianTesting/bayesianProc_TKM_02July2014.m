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
datafile= '/Users/tyler/Dropbox/Coding/Nick_Coding/Data/ducttape_n512_nE32.tnt';


T2guess = [0.02]; %s
nEchoes = 32;       
tEcho   = 150e-6;   % echo time (s)
SNR     = 50;      % approximate SNR
nTrials = 15;        % times to iterate noise
rangeMin = 1e-5;    % minimum examination range (s)
rangeMax = 2e-3;    % maximum examination range (s)
D = 0;         % diffusion constant m^2/s
G = 0;         % field gradient T/m

% Read and set up data file
[aq,spec,~] = readTecmag(datafile);
spec2 = reshape(spec,length(spec)/nEchoes,nEchoes);
s = sum(real(spec2),1)';
s = s./max(s);

% Define prior (in terms of R1 and R2)
if length(T2guess) == 1
    opts.priorMean = [1; 1/T2guess; pi/2; 0];
elseif length(T2guess) == 2
    opts.priorMean = [0.5; 0.5; 1/T2guess(1); 1/T2guess(2); pi/2; 0];
end
opts.priorCov = diag([1000^2*ones(length(T2guess),1); 1000000^2*ones(length(T2guess),1); 100^2*ones(2,1)]);
opts.tau = tEcho;
opts.sigma2 = (1./SNR)^2;
opts.D = D;
opts.G = G;


% Define correction schedule
gammas = logspace(-6,1,30);
gammas = gammas ./ sum(gammas);
opts.correctionSchedule = gammas;

% Allocate output arrays
T2Hat = zeros(length(T2guess),nTrials);
wHat = zeros(length(T2guess),nTrials);

% Repeat for different noise realisations
parfor iTrial=1:nTrials
    
    % Add noise to signal
    y = s + (1./SNR)*randn(size(s));

    % Estimate distribution components
    [T2Hat(:,iTrial), wHat(:,iTrial)] = bayesian_diff_estimate(y,opts);

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