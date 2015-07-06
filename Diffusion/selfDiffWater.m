function [D,Derr,secondary] = selfDiffWater(T_C)

T = T_C + 273.15;

D0 = 1.635e-8; %m^2 s^-1
D0_err = 2.242e-11;

Ts = 215.05; %K
Ts_err = 1.20; 

gamma = 2.063; % unitless
gamma_err = 0.051; 


D = D0*((T/Ts)-1)^gamma;

dD_dD0 = (T/Ts - 1)^gamma;

dD_dTs = -(D0*T*gamma*(T/Ts - 1)^(gamma - 1))/Ts^2;

dD_dgamma = D0*log10(T/Ts - 1)*(T/Ts - 1)^gamma;

Derr = sqrt((dD_dD0)^2*D0_err^2 + (dD_dTs)^2*Ts_err^2+ (dD_dgamma)^2*gamma_err^2);

C1 = [5.8765; 5.7064; 5.3193; 5.6991; 5.6286; 8.4847];
C2 = -[1.6477; 1.6754; 1.6483; 1.7927; 1.8726; 2.9018];
secondary = exp(C1+C2*(1000/T));