clear
clc
close all

cd('/Users/tyler/Desktop/');
load('BSA Data 2.mat');


mm = 66463; %g/mol
samplesPerSec = 5;
rate = 0.8; %deg C/min
mass = 0.259e-3; %g


degPerPt = rate/60/samplesPerSec;
time_minutes = data(:,1);
time_seconds = time_minutes*60;
T = data(:,2);
T_K = T + 273.15;
power = data(:,3); %mW
energy = power/samplesPerSec/1e6; %kJ


meanlims = [1 1700;
            10790 11080];

intlims = [round(mean(meanlims(1,:))) round(mean(meanlims(2,:)))];
delH_m = energy/mass*mm; %kJ mol-1
Cp_m = delH_m./mean(diff(T_K));

figure(4)
hold on
plotyy(T,delH_m,T,Cp_m)



%%



delH_m_start = mean(delH_m(meanlims(1,1):meanlims(1,2)));
delH_m_end = mean(delH_m(meanlims(2,1):meanlims(2,2)));
[Cp_m_max,I] = min(delH_m);
T_maxTrans = T(I);
A = delH_m_end - delH_m_start;
y0 = delH_m_start;

baseline = y0 + A.*(1./(1+exp(0.33*(-T+T_maxTrans))));
subData = -(delH_m - baseline);
delH = sum(subData(intlims(1):intlims(2)))*degPerPt; %deltaH of unfolding in kJ/mol
delHvh = cumsum(subData(intlims(1):intlims(2)))*degPerPt/delH;
K = delHvh./(1-delHvh);

figure(7)
hold on
plot(1./T_K(intlims(1)+2958:intlims(2)),log(K(2959:end)));
plot(T_K(intlims(1)+2959:intlims(2)),8.3145*1e3*abs(diff(1./T_K(intlims(1)+2958:intlims(2))))./(diff(log(K(2959:end)))));

figure(1)
hold on
% plot(t,Cp);
plot(T,delH_m);
plot(T,baseline)
legend('data','sigmoidal baseline')
xlabel('temp (\circC)')
ylabel('enthalpy (kJ/mol)')

figure(2)
hold on
plot(T,subData)
xlim([T(intlims(1)) T(intlims(2))])

C_T = Cp_m./(T_K);

figure(3)
plot(T,C_T)
xlabel('temperature (\circC)')
ylabel('C_{p,m}/T (kJ K^{-1} mol^{-1})')

C_T_start = mean(C_T(meanlims(1,1):meanlims(1,2)));
C_T_end = mean(C_T(meanlims(2,1):meanlims(2,2)));
[C_T_max,I] = min(C_T);
C_T_maxTrans = T(I);
A = C_T_end - C_T_start;
y0 = C_T_start;

baseline = y0 + A.*(1./(1+exp(0.33*(-T+T_maxTrans))));
subC_T = -(C_T - baseline);
delS = sum(subC_T(intlims(1):intlims(2)))*degPerPt*1000/463; %J K-1 (mol residues)-1


T_delG_zero = delH/delS-273.15;

figure(6)
plot(T,subC_T)
xlim([T(intlims(1)) T(intlims(2))])
