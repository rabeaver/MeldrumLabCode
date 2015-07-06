
%% T1T2

cd('Z:\Data\VJL\PM5\1P_2_T2D\1');
dataIn = load('Jamestown_13Nov2014_T1T2data_inverted.out');

T2T1lim = [-3 -1];
T1lim = [-6 1];

T2T1axis = logspace(T2T1lim(1),T2T1lim(2), size(dataIn,1));
T1axis = logspace(T1lim(1),T1lim(2), size(dataIn,2));

figure(1)
surf(T2T1axis,T1axis,dataIn)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('T_2 times (s)')
ylabel('T_1 times (s)')
shading flat


%% T2D
cd('Z:\Data\VJL\PM5\20B_2_T2D\2');
dataInD = load('20B_2_T2D_2.out');

T2Dlim = [-5 0];
Dlim = [-12 -6];

T2Daxis = logspace(T2Dlim(1),T2Dlim(2), size(dataInD,1));
Daxis = logspace(Dlim(1),Dlim(2), size(dataInD,2));

figure(2)
surf(T2Daxis,Daxis,dataInD)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('T_2 times (s)')
ylabel('D (m^2 s^{-1})')
shading flat

%% Combined

figure(3)
subplot(1,2,1)
surf(T2T1axis,T1axis,dataIn)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('T_2 times (s)')
ylabel('T_1 times (s)')
shading flat
subplot(1,2,2)
surf(T2Daxis,Daxis,dataInD)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('T_2 times (s)')
ylabel('D (m^2 s^{-1})')
shading flat