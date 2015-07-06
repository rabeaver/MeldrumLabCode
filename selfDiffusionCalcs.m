clear
close all
clc
%%
T = 20.6 %C

[D, Derr, secondary] = selfDiffWater(T);
% "secondary" has D (x10-9 m2 s-1) for six liquids: cyclohexane, dixane,
% dodecane,DMSO, tetradecane, and pentanol

% temps = 1:1:100;
% 
% for i = 1:length(temps)
%     [D(i),Derr(i)] = selfDiffWater(temps(i));
% end
% 
% figure
% hold on
% plot(temps,D*1e9)
% plot(temps,(D+Derr)*1e9,'--b')
% plot(temps,(D-Derr)*1e9,'--b')
% xlabel('T (C)')
% ylabel('self-diffusion coefficient (x 10^{-9})')

%%
[~,~,secondary] = selfDiffWater(25);
% "secondary" has D (x10-9 m2 s-1) for six liquids: cyclohexane, dixane,
% dodecane,DMSO, tetradecane, and pentanol