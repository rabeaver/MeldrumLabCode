clear
clc
close all

%%

T2_L = 110; %ms, T2 for the ligand only (Nap or Ibu)
T2_B1 = 5; %ms, T2 for the longer relaxation time for BSA
T2_B2 = 0.8; %ms, T2 for the shorter relaxation time for BSA
f_L = 0.5; %fraction of the signal produced by the ligand. (This will be the 0, 0.1, 0.3, etc...)

KD1 = 3.27e-8; %M, KD for the first binding site
KD2 = 5e-5; %M, KD for the second binding site
k_on = 1e10; %s-1, the difusion-limited association of ligand with the host. Probably between 1e7 and 1e10 s-1.

tE = 60e-6; %s, echo time
nE = 1024; %number of echoes
%%

R = diag(1000./[T2_L, T2_B1, T2_B2]);
M0 = [f_L, (1-f_L)/2, (1-f_L)/2]; %initial signal from three sites (L, BSA1, BSA2)

K = [0,         -k_on,      -k_on;
    -KD1*k_on,  0,          0;
    -KD2*k_on,  0,          0];

RK = R + K;
[Q,D] = eig(RK,'nobalance');
RK_arg = -inv(Q)*RK*Q;

%%

echoVec = linspace(tE, nE*tE, nE);
s = zeros(length(M0),length(echoVec));

size(Q)
size(exp(-RK_arg*echoVec(ii)))
size(inv(Q))
size(M0')

for jj = 1:length(M0);
    for ii = 1:length(echoVec)
        s(:,ii) = Q*exp(-RK_arg*echoVec(ii))/Q*(M0');
    end
end

sSum = sum(s,1);
plot(echoVec,sSum)

%%

echoVec = linspace(tE, nE*tE, nE);
s = zeros(length(M0),length(echoVec));

for jj = 1:length(M0);
    for ii = 1:length(echoVec)
        s(jj,ii) = exp(-R(jj,jj)*echoVec(ii))*(M0(jj));
    end
end

sSum = sum(s,1);
plot(echoVec,sSum)

%%
