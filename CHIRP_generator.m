% CHIRP phase table generator for Tecmag
% Adapted from Tecmag's CHIRP script

clear
clc
close all

%%%%%%% User-defined parameters %%%%%%%
% tic
dt = 100e-9 ; %time per point in waveform (s) [Scout limit is 20ns]
% N = 500000; %points to define the waveform
tau = 0.0012; %pulse length (s)
sliceheight = 0.350; %mm
G = 6.59; %T m-1, B0 field gradient
offset = 0; %mm, frequency offset (if applicable)
amplitude = 20; %pwr, for Tecmag
% NOTE: positive offset moves to the left in the FT spectrum (negative
% position)
WURSTshape = 0;
LINEARshape = 1;
linearPct = 0.05; % percent of the front end and back end of the pulse that will be linearly ramped


%%%%%%% END User-defined parameters %%%%%%%


N = round(tau/dt); %number of points per pulse waveform
% dt = (tau/N); %time per point in waveform (s) [Scout limit is 20 ns]
gamma = 42.576; %MHz T-1
SW = sliceheight*G*gamma*1000; %Hz
offsetHz = offset*1000*G*gamma; %Hz
R = SW / tau;

phase = zeros(N,1);
ind = 1:N;
f = (ind-N/2)*(SW/N);
inc = 2*pi*f*dt;
cinc = cumsum(inc);
o_inc = 2*pi*offsetHz*dt;
phase = phase + cinc' - o_inc;
phase = mod((phase*360/2/pi),360);
% toc
dlmwrite('CHIRP_Phase.dat',phase);
% toc

figure(1)
t = linspace(0,tau,N);
plot(t,phase,'-k')
xlabel('time [s]')
ylabel('phase [deg]')
ylim([0 360])


tamp = linspace(0,tau,N);

% WURST shaping
if WURSTshape == 1
    w_amp = amplitude*(1-(cos(pi*tamp/tau)).^40);
    dlmwrite('CHIRP_wAmp.dat',w_amp');
    figure(2)
    plot(tamp,w_amp,'-k')
    xlabel('time [s]')
    ylabel('amplitude [dB]')
    ylim([0 amplitude+1])

% Linear shaping
elseif LINEARshape == 1
    lshapeLeft = linspace(0,1,round(length(tamp)*linearPct));
    lshapeRight = linspace(1,0,round(length(tamp)*linearPct));
    lshape = [lshapeLeft,ones(1,length(tamp)-2*length(lshapeLeft)),lshapeRight];
    l_amp = amplitude*lshape;
    dlmwrite('CHIRP_lAmp.dat',l_amp');
    figure(2)
    plot(tamp,l_amp,'-k')
    xlabel('time [s]')
    ylabel('amplitude [dB]')
    ylim([0 amplitude+1])

end

