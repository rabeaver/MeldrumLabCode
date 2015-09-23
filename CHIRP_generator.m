% CHIRP phase table generator for Tecmag
% Adapted from Tecmag's CHIRP script

clear
clc
close all

%%%%%%% User-defined parameters %%%%%%%
% tic
<<<<<<< HEAD
dt = 1e-6 ; %time per point in waveform (s) [Scout limit is 20ns]
tau = 0.001; %pulse length (s)
sliceheight = 0.100; %mm
=======
dt = 60e-9 ; %time per point in waveform (s) [Scout limit is 20ns]
% N = 500000; %points to define the waveform
tau = 0.008; %pulse length (s)
sliceheight = 0.350; %mm
>>>>>>> 3786e25e63e2af8502f3d632638e3be7b89bb323
G = 6.59; %T m-1, B0 field gradient
offset = 0; %mm, frequency offset (if applicable)
amplitude = 20; %pwr, for Tecmag
% NOTE: positive offset moves to the left in the FT spectrum (negative
% position)

% frequency ramping for CHIRP
LINramp = 1;
EXPramp = 0; %NOT YET FUNCTIONAL

% shape for edges of amplitude profile
WURSTshape = 0;
LINEARshape = 1;
linearPct = 0.05; % percent of the front end and back end of the pulse that will be linearly ramped


%%%%%%% END User-defined parameters %%%%%%%

N = round(tau/dt); %number of points per pulse waveform
gamma = 42.576; %MHz T-1
SW = sliceheight*G*gamma*1000; %Hz
offsetHz = offset*1000*G*gamma; %Hz
t = dt:dt:tau; %time axis for CHIRP pulse
f0 = -SW/2; %initial frequency
    
    
% CALCULATION OF PHASE %

% Linear phase ramping
if LINramp == 1;
    R = SW / tau; %Linear ramp rate
    f = f0 + R*t; %frequency
    inc = 2*pi*f*dt;
    cinc = cumsum(inc);
    o_inc = 2*pi*offsetHz*dt;
    phase = cinc' - o_inc;
    phase = mod((phase*360/2/pi),360);  

% Exponential phase ramping
elseif EXPramp == 1;
    ft = logspace(0,log10(SW),N);
    f = ft - 1 + f0;
    o_inc = 2*pi*offsetHz*dt;
    phase = f - o_inc;
    phase = mod((phase*360/2/pi),360);  
end

%dlmwrite('CHIRP_Phase.dat',phase);


figure(1)
t = linspace(0,tau,N);
plot(t,phase,'-k')
xlabel('time [s]')
ylabel('phase [deg]')
ylim([0 360])


tamp = linspace(0,tau,N);

% CALCULATION OF AMPLITUDE %

% WURST shaping
if WURSTshape == 1
    w_amp = amplitude*(1-(cos(pi*tamp/tau)).^40);
    %dlmwrite('CHIRP_wAmp.dat',w_amp');
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
    %dlmwrite('CHIRP_lAmp.dat',l_amp');
    figure(2)
    plot(tamp,l_amp,'-k')
    xlabel('time [s]')
    ylabel('amplitude [dB]')
    ylim([0 amplitude+1])

end

