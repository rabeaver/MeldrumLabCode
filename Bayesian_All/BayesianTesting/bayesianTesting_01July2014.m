close all
clear
clc

%%
[ap,spec] = readTecmag('/Users/tyler/Dropbox/Coding/Nick_Coding/NAU/glycerol_liquid_n32_19June2014.tnt');

nPts = 32;
nEchoes = length(spec)/nPts;
spec2 = abs(reshape(spec,nPts,nEchoes));
intData = sum(spec2,1);
tau = 150e-6; %echo time in s
t = (1:nEchoes)*tau;
G = 281; %MHz/m
gamma = 42.5; %MHz/T
gammarad = 2.68e8; %rad s-1 T-1

kg = gammarad * tau *G/gamma  ; %rad m-1
D = 1.0e-11; %m^2/s

plot(t,intData)
%%
opts.sigma2 = 0.001;
opts.tau = tau;
gammas = logspace(-6,1,30);
gammas = gammas ./ sum(gammas);
opts.correctionSchedule = gammas;
opts.priorMean = [1; 1/20e-3; pi/2; 0];

nEchoes = length(intData);
N = 2*nEchoes;                % Number of measurements (real & imag)
nParams = length(opts.priorMean);    % Number of parameters to estimate
nComponents = (nParams-2)/2;  % Number of components

opts.priorCov = diag([50^2*ones(nComponents,1); 100000^2*ones(nComponents,1); 100^2*ones(2,1)]);

nEchoes = length(intData);
N = 2*nEchoes;                % Number of measurements (real & imag)
nParams = length(opts.priorMean);    % Number of parameters to estimate
nComponents = (nParams-2)/2;  % Number of components

%%
Findex = 0;	% index of states, 0...N-1
  bvalZ = ((Findex)*kg).^2*tau;		% diffusion  for Z states, assumes that
					% the Z-state has to be refocused, so
					% this models "time between gradients"
	
	% For F states, the following models the additional diffusion time
	% (Findex) and the fact that the state will change if the gradient is
	% on (0.5*Gon), then the additional diffusion *during* the gradient, 
	% ... Gon*kg^2/12 term.

  bvalp = ((( Findex+.5)*kg).^2+kg^2/12)*tau;	% for F+ states
  bvalm = (((-Findex+.5)*kg).^2+kg^2/12)*tau;	% for F- states
  
%%  

% [H, T2, w, flip, phase, covariance]= bayesian_estimate(intData',opts);
% [ H ] = epg( nEchoes, tau, 0, 1/20e-3, pi/2, bvalp, bvalm, bvalZ, D );

H = zeros(nEchoes,1); %one T2 component
alpha = pi/2;
R1 = 0;

T0 = [cos(alpha/2)^2, sin(alpha/2)^2, sin(alpha); ...
    sin(alpha/2)^2, cos(alpha/2)^2, -sin(alpha); ...
    -0.5*sin(alpha), 0.5*sin(alpha), cos(alpha)];

TArray = cell(1,nEchoes);
TArray(:) = {sparse(T0)};
T = blkdiag(1,TArray{:});

% Selection matrix to move all traverse states up one coherence level
%
S = sparse(3*nEchoes+1,3*nEchoes+1);
S(1,3)=1;
S(2,1)=1;
S(3,6)=1;
S(4,4)=1;
for o=2:nEchoes
    offset1=( (o-1) - 1)*3 + 2;
    offset2=( (o+1) - 1)*3 + 3;
    if offset1<=(3*nEchoes+1)
    S(3*o-1,offset1)=1;  % F_k <- F_{k-1}
    end
    if offset2<=(3*nEchoes+1)
    S(3*o,offset2)=1;  % F_-k <- F_{-k-1}
    end
    S(3*o+1,3*o+1)=1;  % Z_order
end

iRate = 1;
    % Relaxation matrix
    R2=1/20e-3;
    R0 = diag([exp(-tau*R2),exp(-tau*R2),exp(-tau*R1)]);

    RArray = cell(1,nEchoes);
    RArray(:) = {sparse(R0)};
    R = blkdiag(exp(-tau*R2),RArray{:});
    
     % Diffusion Matrix
    D0 = diag([exp(-bvalp*D),exp(-bvalm*D),exp(-bvalZ*D)]);
    
    DArray = cell(1,nEchoes);
    DArray(:) = {sparse(D0)};
    D1 = blkdiag((bvalp*D),DArray{:});

    % Precession and relaxation matrix
%     P = (D1*R*S);
    P = (R*S);

    % Matrix representing the inter-echo duration
    E = (P*T*P);
    
    % Recursively apply matrix to get echo amplitudes
    %
    x = zeros(size(R,1),1);
    x(1)=1;
    for iEcho=1:nEchoes
        x=E*x;
        H(iEcho,iRate) = x(1);
    end




%%
figure(2)
hold on
plot(t,abs(H)./max(abs(H)),'-r')
plot(t,intData./max(intData),'-b')