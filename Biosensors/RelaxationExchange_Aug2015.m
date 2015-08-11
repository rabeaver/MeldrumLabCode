clear
clc
close all

%%

T2_L = 110; %ms, T2 for the ligand only (Nap or Ibu)
T2_B1 = 5; %ms, T2 for the longer relaxation time for BSA
T2_B2 = 0.8; %ms, T2 for the shorter relaxation time for BSA
f_L = 1; %fraction of the signal produced by the ligand. (This will be the 0, 0.1, 0.3, etc...)

KD1 = 3.27e-8; %M, KD for the first binding site
KD2 = 5e-5; %M, KD for the second binding site
k_on = 1e7; %s-1, the difusion-limited association of ligand with the host. Probably between 1e7 and 1e10 s-1.

tE = 600e-6; %s, echo time
nE = 128; %number of echoes
%%
R = diag(1000./[T2_L,T2_B1]); %relaxation rate matrix
M0 = ([ f_L; 1-f_L]); %initial magnetization (really, initial signal from each component)


% K = [KD1*k_on,     -k_on;
%      -KD1*k_on,     k_on];

K = [1000*(1-f_L)/f_L, -1000; %exchange matrix, need to be sure to balance magnetization
     -1000*(1-f_L)/f_L, 1000];

% this part determines the eigenvals, eigenvectors, and rotation matrices
% for the R+K matrix
RK = R + K; 
[Q,D] = eig(RK);
RK_arg = -Q\RK*Q;

%%
echoVec = linspace(tE, nE*tE, nE); %make vector for time axis of echoes
s = zeros(length(M0),length(echoVec)); %holder for signal


    for ii = 1:length(echoVec)
        s(:,ii) = Q*exp(RK_arg*echoVec(ii))/Q*(M0)+M0; %main calculation
    end

sSum = sum(s,1);
sSum = sSum./sSum(1);
plot(echoVec,sSum)
