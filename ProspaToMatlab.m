
%% T1T2

cd('C:\users\vjlee\desktop\');
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
cd('Z:\Data\JNK\PM5\Paint_Nick\1July2016_to_4July2016\M212_2014_TradSSET2_Overnight_4July2016\1\');
dataInD = load('ProspaILT_correctaxes.out');

T2Dlim = [-5 0];
Dlim = [-12 -7];

T2Daxis = logspace(T2Dlim(1),T2Dlim(2), size(dataInD,1));
Daxis = logspace(Dlim(1),Dlim(2), size(dataInD,2));

figure(2)
surf(T2Daxis,Daxis,dataInD)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('T_2 (s)')
ylabel('D (m^2 s^{-1})')
shading flat
view(0,90)

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