datadir = '/Users/tyler/Desktop/TOSEND_Genpurp_DELM_15Aug/1/';
datafile = 'dataRe.dat';
paramsfile = 'acqu.par';

params.acqTime = readpar_Kea(strcat(datadir,paramsfile),'acqTime');
params.bandwidth = readpar_Kea(strcat(datadir,paramsfile),'bandwidth');
params.nScans = readpar_Kea(strcat(datadir,paramsfile),'nrScans');
params.rxPhase = readpar_Kea(strcat(datadir,paramsfile),'rxPhase');
params.rxGain = readpar_Kea(strcat(datadir,paramsfile),'rxGain');
params.nPts = readpar_Kea(strcat(datadir,paramsfile),'nrPnts');
params.repTime = readpar_Kea(strcat(datadir,paramsfile),'repTime');
params.b1Freq = readpar_Kea(strcat(datadir,paramsfile),'b1Freq');
params.nEchoes = readpar_Kea(strcat(datadir,paramsfile),'nrEchoes');
params.tE = readpar_Kea(strcat(datadir,paramsfile),'echoTime');
