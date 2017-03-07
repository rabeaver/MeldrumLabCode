% CHIRP phase table generator for Tecmag
% Adapted from Tecmag's CHIRP script

clear
clc
close all

%%%%%%% User-defined parameters %%%%%%%
% tic


dt = 1e-6;             % time per point in waveform (s) [Scout limit is 20ns]
tau = 300e-3;            % pulse length (s)
sliceheight = 0.250;     % mm


G = 6.59;               %T m-1, B0 field gradient [PM25 = 6.59]
                        %                         [PM5 = 23.87]
offset = 0;             %mm, frequency offset (if applicable)
amplitude = 5;         %pwr, for Tecmag
% NOTE: positive offset moves to the left in the FT spectrum (negative
% position)

% frequency ramping for CHIRP
LINramp = 1;
EXPramp = 0; 

% shape for edges of amplitude profile
WURSTshape = 0;
LINEARshape = 1;
linearPct = 0.01;       % percent of the front end and back end of the pulse
                        % that will be linearly ramped


%%%%%%% END User-defined parameters %%%%%%%

N = round(tau/dt);              %number of points per pulse waveform
gamma = 42.576;                 %MHz T-1
SW = sliceheight*G*gamma*1000;  %Hz
% SW_s = SW*2*pi;               %rad s-1
offsetHz = offset*1000*G*gamma; %Hz
f0 = -SW/2;                     %initial frequency (Hz)
% f0_s = -SW_s/2;               %initial frequency (s-1)
    
    
% CALCULATION OF PHASE %

% Linear phase ramping
if LINramp == 1;
    
%analytical expressions *POSSIBLY INCORRECT WITH Hz/rad/deg UNITS*
%     f = f0 + SW*t/tau;                          %frequency (Hz, or cycles s-1)
%     f_s = f0_s + SW_s*t/tau;                    %frequency (rad s-1)
%     phase = f0*t+SW/2/tau*t.^2;                 %phase (cycles, or Hz*s)
%     phase_s = f0_s*t + SW_s/2/tau*t.^2;         %phase (rad)

% numerical expressions
    % WE FLIPPED F SHOULD BE LINSPACE(-SW/2,SW/2,N)!!!!!!!!!!!!!!!!!!!
        % WE FLIPPED F SHOULD BE LINSPACE(-SW/2,SW/2,N)!!!!!!!!!!!!!!!!!!!
            % WE FLIPPED F SHOULD BE LINSPACE(-SW/2,SW/2,N)!!!!!!!!!!!!!!!!!!!
                % WE FLIPPED F SHOULD BE LINSPACE(-SW/2,SW/2,N)!!!!!!!!!!!!!!!!!!!
    f = linspace(-SW/2,SW/2,N);       % frequency (Hz, or cycles s-1)
    % WE FLIPPED F SHOULD BE LINSPACE(-SW/2,SW/2,N)!!!!!!!!!!!!!!!!!!!
        % WE FLIPPED F SHOULD BE LINSPACE(-SW/2,SW/2,N)!!!!!!!!!!!!!!!!!!!
            % WE FLIPPED F SHOULD BE LINSPACE(-SW/2,SW/2,N)!!!!!!!!!!!!!!!!!!!
                % WE FLIPPED F SHOULD BE LINSPACE(-SW/2,SW/2,N)!!!!!!!!!!!!!!!!!!!
    f_s = f*360;                      % frequency (deg s-1)
    phase_s = cumsum(f_s)*dt;         % phase (deg) 
    phase_s360 = mod(phase_s,360);    % phase (deg, mod 360)

    
    
% Exponential phase ramping
elseif EXPramp == 1;
    
%analytical expression for exponential growth *POSSIBLY INCORRECT WITH Hz/rad/deg UNITS*
%     f = SW.^(t/(tau*(1-linearPct))) - 1 + f0;
%     phase = (tau*(1-linearPct)/log(SW))*(SW.^(t/(tau*(1-linearPct)))-1)-t*(1-f0);
    
% numerical log-spaced points
    f = logspace(0,log10(SW),N)-1-SW/2;       % frequency (Hz, or cycles s-1)
    f_s = f*360;                              % frequency (deg s-1)
    phase = cumsum(f)*dt;                     % phase (cycles, or Hz*s)
    phase_s = cumsum(f_s)*dt;                 % phase (deg)
    phase_s360 = mod(phase_s,360);            % phase (deg, mod 360)
end

dlmwrite('CHIRP_Phase.dat',phase_s360');

t = linspace(0,tau,N);
figure(1)
subplot(1,2,1)
plot(t,f_s,'-k')
xlabel('time [s]')
ylabel('frequency [s-1]')
subplot(1,2,2)
plot(t,phase_s360,'-k')
xlabel('time [s]')
ylabel('phase [deg]')
ylim([0 360])


tamp = linspace(0,tau,N);

% CALCULATION OF AMPLITUDE %
amp = amplitude; 

% WURST shaping
if WURSTshape == 1
    amp = amplitude*(1-(cos(pi*tamp/tau)).^40);
    
% Linear shaping
elseif LINEARshape == 1
    lshapeLeft = linspace(0,1,round(length(tamp)*linearPct));
    lshapeRight = linspace(1,0,round(length(tamp)*linearPct));
    lshape = [lshapeLeft,ones(1,length(tamp)-2*length(lshapeLeft)),lshapeRight];
    amp = amplitude*lshape; 
end


dlmwrite('CHIRP_Amp.dat',amp');
    figure(2)
    plot(tamp,amp,'-k')
    xlabel('time [s]')
    ylabel('amplitude [pwr]')
    ylim([0 amplitude+1])

        
    figure(3)
    plot(tamp,amp.*phase_s360/360,'-b')
    xlabel('time [s]')
    ylabel('phase * amplitude [dB]')
    ylim([0 amplitude+1])