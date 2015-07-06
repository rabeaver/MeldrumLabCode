clear
cd('Z:\ISC1026\Data\JYU\BSA_Nap\T2_D_BSA_Only\2');
signal = load('dataRe.dat');
data = load('data.dat');
tE = readpar_Kea('acqu.par','echoTime')/1000; % get echo time value from acqu.par file
NE = readpar_Kea('acqu.par','nrEchoes'); % get echo time value from acqu.par file

gamma = 267.513e6; %T-1 s-1
G = 23.6/100; %T cm-1
DELTA = 1e-3; %s

t1 = tE:tE:NE*tE;
t1 = t1/1000;
t2 = data(:,1)/1000;
q = (gamma*G*t2).^2.*(DELTA-(2/3)*t2);

dlmwrite('signal.out',signal);
dlmwrite('t_h.out',t1);
dlmwrite('q_v.out',q);

% Z:\ISC1026\Data\JYU\BSA_Nap\T2_D_NPna_Only\1