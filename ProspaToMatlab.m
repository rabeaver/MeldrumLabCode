
%% T1T2

cd('Z:\ISC1026\Data\TKM\Mortar\T2D');
dataIn = load('OneToOneWater_T2D_17Nov2014_Inverted.out');

T2T1lim = [-3 -1];
T1lim = [-3 -1];

T2T1axis = logspace(T2T1lim(1),T2T1lim(2), size(dataIn,1));
T1axis = logspace(T1lim(1),T1lim(2), size(dataIn,2));

% figure(1)
% surf(T2T1axis,T1axis,dataIn)
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% xlabel('T_2 times (s)')
% ylabel('T_1 times (s)')
% shading flat


%% T2D
cd('Z:\ISC1026\Data\TKM\Mortar\T2D\');
dataInD = load('OneToOne_T2D_14Nov2014_Inverted.out');

T2Dlim = [-5 0];
Dlim = [-9 -7];

T2Daxis = logspace(T2Dlim(1),T2Dlim(2), size(dataInD,1));
Daxis = logspace(Dlim(1),Dlim(2), size(dataInD,2));

figure(2)
pcolor(T2Daxis,Daxis,dataInD)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('T_2 times (s)')
ylabel('D (m^2 s^{-1})')
shading flat
% view(0,90)

%% Combined

figure(3)
subplot(1,2,1)
surf(T1axis,T2T1axis,dataIn')
% line([min(T1axis) min(T1axis)],[max(T1axis) max(T1axis)]);
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('T_1 times (s)')
ylabel('T_2 times (s)')
shading flat
view(0,90)
subplot(1,2,2)
surf(Daxis,T2Daxis,dataInD')
set(gca,'XScale','log')
set(gca,'YScale','log')
ylabel('T_2 times (s)')
xlabel('D (m^2 s^{-1})')
shading flat
view(0,90)