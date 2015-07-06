close all
clear
clc

%% 

% ----- User Defined Parameters
filename = '/Users/tyler/Dropbox/Data/T1Coding/Jsoakedwater_trep8_T1_n21_nE32_s4_3July2014.tnt';
nEchoes = 32;   %number of Echoes
zPts = 5;   %number of zero data points at end of acq window
T1lin = 10; % 1 if linearly spaced tau pts, 0 if logarithmically spaced
tauEst = 1.6e6; %estimated T1 in us
pw = 6; %pulse length in us
% ----- End user-defined Parameters




[ap, spec, spec2] = readTecmag(filename);
nPts = ap.td(1)/nEchoes;    %number of complex points


% set up time points
n = 1:ap.td(2);
if T1lin == 1;
    tau = (0.05*tauEst) + ((n-1).*(4.95*tauEst)./(ap.td(2)-1)) - pw;
    tau = tau';
else
    tau = -tauEst * log(1- ((1-exp(-5))/ap.td(2)).*n)-pw;
    tau = tau';
end


% sums and reshapes the data
spec3 = reshape(spec2',nPts,nEchoes,ap.td(2));
spec3 = sum(spec3(1:nPts-zPts,:,:));
data = reshape(spec3,nEchoes,ap.td(2));
T1data = real(sum(data,1))./max(real(sum(data,1)));

guesses = [0, 1, 1];
beta = nlinfit(tau*1e-6,T1data',@T1_recovery,guesses);
pred = T1_recovery(beta,tau*1e-6);

T1 = beta(3) %T1 value in s

%% plot data
figure(1)
hold on
plot(tau*1e-6,pred,'-r')
scatter(tau*1e-6,T1data,'b')
ylim([0 1.1])
xlabel('recovery time [s]')
ylabel('signal [arb]')
legend('fit','data')