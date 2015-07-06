function [contrast]=contrastcalc(tauoff,tauon,time)
contrast = (exp((-time)./tauoff) - exp((-time)./tauon))./(exp((-time)./tauoff));