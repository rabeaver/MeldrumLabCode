clear
clc
close all

%%
[aq,spec,spec2] = readTecmag('C:\Users\tkmeldrum\Dropbox\Data\NAU\glycerol_liquid_n32_30June2014_T1.tnt');
intData = sum(abs(spec2),2);

tEst = 17000; %us
nPts = 11;
pw = 6e-6;
logYN = 1;

if logYN == 0
    t = (0.05*tEst)+((1:nPts)' - 1).*(4.95*tEst/(nPts - 1)) - pw;
elseif logYN == 1
    t1 =  1 - (((1-exp(-5))/(nPts)) .* (1:nPts)');
    t = -tEst .* log(t1) - pw;
end


