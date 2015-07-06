clear
clc
close all

T2{1} = [30.2; 2.4];
D{1} = [1.5e-10; 1.04e-10];
A{1} = [0.37; 0.63];

T2{2} = [28.2; 3.16; 0.316];
D{2} = [1.51e-10; 8.71e-11; 8.32e-11];
A{2} = [0.32; 0.60; 0.085];

T2{3} = [24.2; 2.65; 0.0912];
D{3} = [1.18e-10; 1.12e-10; 2.13e-10];
A{3} = [0.53; 0.38; 0.084];

T2{4} = [25.1; 1.12];
D{4} = [1.15e-10; 8.32e-11];
A{4} = [0.69; 0.29];

T2{5} = [48.6; 2.51];
D{5} = [1.08e-10; 1.09e-10];
A{5} = [0.82; 0.18];

T2{6} = [63.1; 0.562];
D{6} = [1e-10; 6.92e-11];
A{6} = [0.86; 0.039];

T2{7} = [70.8];
D{7} = [9.55e-11];
A{7} = [0.67];

BSAT2{1} = [1.19; 7.08];
BSAD{1} = [8.91e-11; 8.91e-11];
BSAA{1} = [0.5; 0.5];

BSAT2{2} = [4.17];
BSAD{2} = [9.12e-11];
BSAA{2} = [0.79];

BSAT2{3} = [4.17];
BSAD{3} = [8.63e-11];
BSAA{3} = [1];

NapT2{1} = [89.8; 7.59];
NapD{1} = [1.08e-10; 1.22e-10];
NapA{1} = [0.69; 0.67];

NapT2{2} = [40.4];
NapD{2} = [1.25e-10];
NapA{2} = [0.48];

NapT2{3} = [35.3];
NapD{3} = [5.97e-10];
NapA{3} = [0.24];


%%
n = 7;

% close all

if n == 1;
    titleText = 'BSA 10:0 Naproxen';
elseif n == 2;
    titleText = 'BSA 9:1 Naproxen';
elseif n == 3;
    titleText = 'BSA 7:3 Naproxen';
elseif n == 4;
    titleText = 'BSA 5:5 Naproxen';
elseif n == 5;
    titleText = 'BSA 3:7 Naproxen';
elseif n == 6;
    titleText = 'BSA 1:9 Naproxen';
elseif n == 7;
    titleText = 'BSA 0:10 Naproxen';
end

figure(n)
hold on
scatter(T2{n},D{n},500*A{n},'filled')
set(gca,'yscale','log','xscale','log')
xlim([0.01 100])
ylim([1e-11 1e-9])
xlabel('T_2 [ms]')
ylabel('D [m^2 s^{-1}]')
title(titleText)


%%
n = 9;

if n == 1;
    titleText = 'BSA 10 [:0 Naproxen]';
    figure(n)
    scatter(T2{1},D{1},500*A{1},'filled')
    set(gca,'yscale','log','xscale','log')
    xlim([0.01 100])
    ylim([1e-11 1e-9])
    xlabel('T_2 [ms]')
    ylabel('D [m^2 s^{-1}]')
    title(titleText)
elseif n == 2;
    titleText = 'BSA 9 [:1 Naproxen]';
    figure(n)
    scatter(BSAT2{1},BSAD{1},500*BSAA{1},'filled')
    set(gca,'yscale','log','xscale','log')
    xlim([0.01 100])
    ylim([1e-11 1e-9])
    xlabel('T_2 [ms]')
    ylabel('D [m^2 s^{-1}]')
    title(titleText)
elseif n == 3;
    titleText = 'BSA 5 [:5 Naproxen]';
    figure(n)
    scatter(BSAT2{2},BSAD{2},500*BSAA{2},'filled')
    set(gca,'yscale','log','xscale','log')
    xlim([0.01 100])
    ylim([1e-11 1e-9])
    xlabel('T_2 [ms]')
    ylabel('D [m^2 s^{-1}]')
    title(titleText)
elseif n == 4;
    titleText = 'BSA 1 [:9 Naproxen]';
    figure(n)
    scatter(BSAT2{3},BSAD{3},500*BSAA{3},'filled')
    set(gca,'yscale','log','xscale','log')
    xlim([0.01 100])
    ylim([1e-11 1e-9])
    xlabel('T_2 [ms]')
    ylabel('D [m^2 s^{-1}]')
    title(titleText)
elseif n == 5;
    titleText = '[BSA 0:] 10 Naproxen';
    figure(n)
    scatter(T2{7},D{7},500*A{7},'filled')
    set(gca,'yscale','log','xscale','log')
    xlim([0.01 100])
    ylim([1e-11 1e-9])
    xlabel('T_2 [ms]')
    ylabel('D [m^2 s^{-1}]')
    title(titleText)
elseif n == 6;
    titleText = '[BSA 1:] 9 Naproxen';
    figure(n)
    scatter(NapT2{1},NapD{1},500*NapA{1},'filled')
    set(gca,'yscale','log','xscale','log')
    xlim([0.01 100])
    ylim([1e-11 1e-9])
    xlabel('T_2 [ms]')
    ylabel('D [m^2 s^{-1}]')
    title(titleText)    
elseif n == 7;
    titleText = '[BSA 5:] 5 Naproxen';
    figure(n)
    scatter(NapT2{2},NapD{2},500*NapA{2},'filled')
    set(gca,'yscale','log','xscale','log')
    xlim([0.01 100])
    ylim([1e-11 1e-9])
    xlabel('T_2 [ms]')
    ylabel('D [m^2 s^{-1}]')
    title(titleText)   
elseif n == 8;
    titleText = '[BSA 9:] 1 Naproxen';
    figure(n)
    scatter(NapT2{3},NapD{3},500*NapA{3},'filled')
    set(gca,'yscale','log','xscale','log')
    xlim([0.01 100])
    ylim([1e-11 1e-9])
    xlabel('T_2 [ms]')
    ylabel('D [m^2 s^{-1}]')
    title(titleText)   
elseif n == 9;
    titleText = '[BSA 10:] 0 Naproxen';
    figure(n)
    scatter(T2{1},D{1},500*A{1},'filled')
    set(gca,'yscale','log','xscale','log')
    xlim([0.01 100])
    ylim([1e-11 1e-9])
    xlabel('T_2 [ms]')
    ylabel('D [m^2 s^{-1}]')
    title(titleText)
end

% figure(1)
% scatter(T2{n},D{n},500*A{n},'filled')
% set(gca,'yscale','log','xscale','log')
% xlim([0.01 100])
% ylim([1e-11 1e-9])
% xlabel('T_2 [ms]')
% ylabel('D [m^2 s^{-1}]')
% title(titleText)
