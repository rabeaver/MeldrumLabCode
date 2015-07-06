clc 
close all
clear 

% cd('Z:\BMF\TdExpParams')% wherever you want the text files saved
cd('Z:\BMF\TdExpParams')
%% Excel TD diff time calculator

% params
tDpts = 12; % for varying total TD
tDmin = 60;%in us,
tDmax = 460;%in us, 

nrDpts = 9; % per TD calc
tEcho = 150; % us
dwellTime = 1; %in us
acqPts = 69; 
pulseLength = 6; %us 
% deltaMin = 7; % us for fixed deltas must be commented out below, otherwise its 0.1 *td for min, 0.4* Td for max
% deltaMax = 31; % us for fixed deltas
nrScans = 512;
tRep = 1;% in s

experimentName = ('Vrestricted_42315_TD');

% calcs
tD = linspace(tDmin,tDmax,tDpts); %in us,
expTime = tDpts*nrDpts*nrScans*tRep/3600 % gives exp time in hours

% deltaSpace = (deltaMax - deltaMin)./(nrDpts-1);
steps = 1:nrDpts;
acqTime = acqPts*dwellTime; % us

for j = 1:tDpts;% goes to the number of tD points
    
deltaMin = 0.1*tD(j); % calculates Deltas based on the tD
deltaMax = 0.4*tD(j); % calculates Deltas based on the tD
deltaSpace = (deltaMax - deltaMin)./(nrDpts-1);
    
for i =steps
    delta(i) = deltaMin+deltaSpace*((steps(i)-1));
    DELTA(i) = tD(j)-2*delta(i);
    d1(i) = delta(i) - pulseLength;
    d2(i) = DELTA(i) - pulseLength;
    d3(i) = delta(i) + tEcho/2 - pulseLength/2;
    d4(i) = tEcho/2 - acqTime/2 -pulseLength/2;
end
nrpts(j).delta = delta;
nrpts(j).DELTA = DELTA;
nrpts(j).d1 = d1;
nrpts(j).d2 = d2;
nrpts(j).d3 = d3;
nrpts(j).d4 = d4;

end

%% concatenating points and printing into a text file for use
delta = nrpts(1).delta';  %' so that the text file will be a large column
DELTA = nrpts(1).DELTA';
d1 = nrpts(1).d1';
d2 = nrpts(1).d2';
d3 = nrpts(1).d3';
d4 = nrpts(1).d4';

for j = 2:tDpts;
    delta = [delta; nrpts(j).delta'];
    DELTA = [DELTA; nrpts(j).DELTA'];
    d1 = [d1; nrpts(j).d1'];
    d2 = [d2; nrpts(j).d2'];
    d3 = [d3; nrpts(j).d3'];
    d4 = [d4; nrpts(j).d4'];
end

save((strcat(experimentName,'delta.txt')),'delta','-mat');
save((strcat(experimentName,'DELTA.txt')),'DELTA','-ascii');
save((strcat(experimentName,'d1.txt')),'d1','-ascii');
save((strcat(experimentName,'d2.txt')),'d2','-ascii');
save((strcat(experimentName,'d3.txt')),'d3','-ascii');
save((strcat(experimentName,'d4.txt')),'d4','-ascii');
