function [alpha90, alpha180] = findalpha(t90, dB90, tCHIRP, swCHIRP, G, CHIRP90pwr, CHIRP180pwr)
% Function: findalpha
% Author: Jared King
% Date: 06/06/2016
%
% This function finds the value of the constants 'alpha' for the CHIRP
% pulse powers
%
% Variables Passed: t90, dB90, tCHIRP, swCHIRP, G, CHIRP90pwr, CHIRP180pwr
% t90           - Calibrated time for hard 90 degree pulse (s)
% dB90          - Power used for hard 90 degree pulse (dB)
% tCHIRP        - Time of CHIRP pulse (s)
% swCHIRP       - Sliceheight of the CHIRP pulse (um)
% G             - The Gradient of the magnet in use (T/m)
% CHIRP90pwr    - Calibrated power of 90 degree CHIRP pulse (dB)
% CHIRP180pwr   - Calibrated power of 180 degree CHIRP pulse (dB)
% Note: Passing a number greater than 0 for CHIRP90pwr or CHIRP180pwr will
% bypass the calculation for that alpha
% e.g. passing findalpha(~,~,~,~,~,1, -60) will find alpha for the 180
% degree CHIRP pulse but not the 90 degree CHIRP pulse
%
%
% Variables Returned: alpha90, alpha180 (units?)
hpPwr = 1e6/t90; %Hz
bwCHIRP

end
