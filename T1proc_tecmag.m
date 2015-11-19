clear
clc
close all

%%
[aq,spec,spec2] = readTecmag('C:\CommonData\CHIRP\T2D\ShellSol_T1IR_BURP_29Oct2015.tnt');
intData = sum(abs(spec2),2);

tEst = 200000; %us
nPts = 11;
pw = 28e-6;
logYN = 0;

if logYN == 0
    t = (0.05*tEst)+((1:nPts)' - 1).*(4.95*tEst/(nPts - 1)) - pw;
elseif logYN == 1
    t1 =  1 - (((1-exp(-5))/(nPts)) .* (1:nPts)');
    t = -tEst .* log(t1) - pw;
end

scatter(t,intData)