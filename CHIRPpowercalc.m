function [pwr90,pwr180] = CHIRPpowercalc(t90,dB90,tCHIRP,bwCHIRP,G,alpha90,alpha180) %t in us, sw in um

hpPwr = 1e6/(4*t90); %Hz
swCHIRP = bwCHIRP*G*42.576;
chirpR = swCHIRP/(tCHIRP*1e-6); %Hz /s
chirp90pwr = alpha90*sqrt(chirpR); %chirp90 power  in Hz
chirp180pwr = alpha180*sqrt(chirpR); %chirp90 power  in Hz

p_c90hp = chirp90pwr/hpPwr; %chirp90 pwr as fraction of hard pulse power
p_c180hp = chirp180pwr/hpPwr; %chirp180 pwr as fraction of hard pulse power

lin_hppwr = round(2^14*10^(dB90/20-1)); %linearized hard pulse power

p_c90 = p_c90hp*lin_hppwr; %linearized chirp 90 power
p_c180 = p_c180hp*lin_hppwr; %linearized chirp 180 power

pwr90 = 20*(log10(p_c90)-14*log10(2)+1); %chirp90 power in dB
pwr180 = 20*(log10(p_c180)-14*log10(2)+1); %chirp180 power in dB
end


