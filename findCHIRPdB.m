function [CHIRPpwr90, CHIRPpwr180] = findCHIRPdB(hp90time, hp90dB, CHIRPtime, swCHIRP, G, alpha90, alpha180)
% Function: findCHIRPdB.m
% Author: Jared King
% Date: 06/09/2016
%
% This function uses the previously determined values of alpha to determine
% pulse powers (dB) of 180 and 90 degree CHIRP pulses
%
% Variables Passed: hp90time, hp90dB, CHIRPtime, swCHIRP, G, alpha90, alpha180
% hp90time      - Calibrated time for hard 90 degree pulse (s)
% hp90dB        - Power used for hard 90 degree pulse (dB)
% CHIRPtime     - Time of CHIRP pulse (s)
% swCHIRP       - Sliceheight of the CHIRP pulse (um)
% G             - The Gradient of the magnet in use (T/m)
% alpha90       - Constant alpha for 90 degree CHIRP pulse
% alpha180      - Constant alpha for 180 degree CHIRP pulse
% Note: Passing a number less than 0 for alpha90 or alpha180 will
% bypass the calculation for that pulse power
% e.g. passing findCHIRPdB(~,~,~,~,~,-1, 5) will find alpha for the 180
% degree CHIRP pulse but not the 90 degree CHIRP pulse
%
%
% Variables Returned: CHIRPpwr90 (dB), CHIRPpwr180 (dB)

hpPwr = 4/hp90time; % Hz
bwCHIRP = swCHIRP * 42.576 * G; %Bandwidth of CHIRP in Hz
CHIRPr = bwCHIRP/(CHIRPtime); % s^-2
dB90 = hp90dB;

if alpha90 ~= -1
    linpwr90 = sqrt(CHIRPr) * alpha90;
    CHIRPpwr90 = (log10(linpwr90/hpPwr) * 20) + dB90;
    
else
    CHIRPpwr90 = NaN;
end

if alpha180 ~= -1
    linpwr180 = sqrt(CHIRPr) * alpha180;
    CHIRPpwr180 = (log10(linpwr180/hpPwr) * 20) + dB90;
    
else
    CHIRPpwr180 = NaN;
end


end