t90 = 12e-6; %s, hard 90 pulse time
p90 = 30; %hard 90 pulse power, linear scale (ecmag)

BW = 56000; %Hz, CHIRP bandwidth
tCHIRP = 200e-6; %s, CHIRP duration

% convert rf pulse time and power to frequency power
t360 = 4*t90;
frf = 1/t360; %Hz

%CHIRP pulse sweep rate
fSweep = BW/tCHIRP; %Hz s-1 (or s-2)

%coefficients from:
% (1) Shrot, Y., & Frydman, L. J. Mag. Reson. 2005, 172(2), 179?90.
% (2) Bhattacharyya, R., & Frydman, L. J. Chem. Phys. 2007, 127(19), 194503.
% These predict the correct conditions for a 90 and 180 CHIRP pulse.
a90 = 0.26;
a180 = 0.8;

%determine rf power for CHIRP pulses in Hz (prf...) and in linear (0-100)
%power (pSetting)
prfCHIRP90 = a90*sqrt(fSweep);
prfCHIRP180 = a180*sqrt(fSweep);

pSetting90 = p90*(prfCHIRP90/frf)
pSetting180 = p90*2*(prfCHIRP180/frf)