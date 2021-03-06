clear
clc
close all

sampleCt = 3;
viewlim = [1e-4 1e-1];
%%
% Glycerol
data = load('/Users/jaredking/Documents/Classes/Chemistry/Research/Summer2016/TradSSET2_M212Samples_23_28June2016/TradSSET2_M212_9months_Overnight_23June2016/1/ILTs/TradSSET2_M212_Overnight_23June2016.out');
lowlim = 1e-4;
hilim = 1e-1;

xaxis = logspace(log10(lowlim), log10(hilim), size(data,1));
yaxis = logspace(log10(lowlim), log10(hilim), size(data,2));

h = figure('Units','centimeters','Position', [0 0 19 28])
subplot(sampleCt,2,1)
grid on
surf(xaxis,yaxis,data./sum(sum(data)))
set(gca,'XScale','log','YScale','log','TickLabelInterpreter', 'latex','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
colormap('gray')
colormap(flipud(colormap))
shading flat
view([0 90])
line([lowlim hilim],[lowlim hilim],[max(max(data)) max(max(data))],'LineWidth',2,'Color','red','LineStyle','--')
xlabel({'$\it{T}\rm{_2 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
ylabel({'$\it{T}\rm{_1 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
title('Glycerol, Traditional','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
colorbar
xlim(viewlim)
ylim(viewlim)

clear('data','xaxis','yaxis')

data = load('/Users/tyler/Desktop/CHIRP_Manuscript/Raw Data/Glycerol/Glycerol_CHIRP_29Sep2015.out');
% lowlim = 1e-4;
% hilim = 1e-2;

xaxis = logspace(log10(lowlim), log10(hilim), size(data,1));
yaxis = logspace(log10(lowlim), log10(hilim), size(data,2));

figure(h)
subplot(sampleCt,2,2)
surf(xaxis,yaxis,data./sum(sum(data)))
% colormap('gray')
colormap(flipud(colormap))
shading flat
view([0 90])
set(gca,'XScale','log','YScale','log','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
line([lowlim hilim],[lowlim hilim],[max(max(data)) max(max(data))],'LineWidth',2,'Color','red','LineStyle','--')
xlabel({'$\it{T}\rm{_2 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
ylabel({'$\it{T}\rm{_1 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
title('Glycerol, CHIRP','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
colorbar
xlim(viewlim)
ylim(viewlim)
%
% 15 mM Gd
data = load('/Users/tyler/Desktop/CHIRP_Manuscript/Raw Data/15mMGd/Inverted_FISTA6000_400x400_15mM_T1IRBURP.out');
lowlim = 1e-4;
hilim = 1e-1;

xaxis = logspace(log10(lowlim), log10(hilim), size(data,1));
yaxis = logspace(log10(lowlim), log10(hilim), size(data,2));

figure(h)
subplot(sampleCt,2,3)
surf(xaxis,yaxis,data./sum(sum(data)))
colormap('gray')
colormap(flipud(colormap))
shading flat
view([0 90])
set(gca,'XScale','log','YScale','log','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
line([lowlim hilim],[lowlim hilim],[max(max(data)) max(max(data))],'LineWidth',2,'Color','red','LineStyle','--')
xlabel({'$\it{T}\rm{_2 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
ylabel({'$\it{T}\rm{_1 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
title('15 mM Gd, noCHIRP','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
xlim(viewlim)
ylim(viewlim)
colorbar

clear('data','xaxis','yaxis')

data = load('/Users/tyler/Desktop/CHIRP_Manuscript/Raw Data/15mMGd/Inverted_FISTA6000_400x400_15mM_CHIRP.out');
% lowlim = 1e-4;
% hilim = 1e-2;

xaxis = logspace(log10(lowlim), log10(hilim), size(data,1));
yaxis = logspace(log10(lowlim), log10(hilim), size(data,2));

figure(h)
subplot(sampleCt,2,4)
surf(xaxis,yaxis,data./sum(sum(data)))
% colormap('gray')
% colormap(flipud(colormap))
shading flat
view([0 90])
set(gca,'XScale','log','YScale','log','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
line([lowlim hilim],[lowlim hilim],[max(max(data)) max(max(data))],'LineWidth',2,'Color','red','LineStyle','--')
xlabel({'$\it{T}\rm{_2 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
ylabel({'$\it{T}\rm{_1 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
title('15 mM Gd, CHIRP','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
xlim(viewlim)
ylim(viewlim)
colorbar
%
% Double
% data = load('/Users/tyler/Dropbox/Manuscripts/CHIRP/ExcelTableData/DoubleSample_noCHIRP.out');
lowlim = 1e-4;
hilim = 1e0;

% xaxis = logspace(log10(lowlim), log10(hilim), size(data,1));
% yaxis = logspace(log10(lowlim), log10(hilim), size(data,2));
% 
% figure(h)
% subplot(sampleCt,2,5)
% surf(xaxis,yaxis,data./sum(sum(data)))
% colormap('gray')
% colormap(flipud(colormap))
% shading flat
% view([0 90])
% set(gca,'XScale','log','YScale','log','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
% line([lowlim hilim],[lowlim hilim],[max(max(data)) max(max(data))],'LineWidth',2,'Color','red','LineStyle','--')
% xlabel({'$\it{T}\rm{_2 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
% ylabel({'$\it{T}\rm{_1 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
% title('Double noCHIRP','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
% xlim(viewlim)
% ylim(viewlim)

clear('data','xaxis','yaxis')

data = load('/Users/tyler/Desktop/CHIRP_Manuscript/Raw Data/Double_15mM_Glycerol/CHIRP_DOUBLE_15mM_Gly_40mspw_sliceheight350um_tD8u_76pts_1024scans_100nsWave_29Sept2015.out');
% lowlim = 1e-4;
% hilim = 1e-2;

xaxis = logspace(log10(lowlim), log10(hilim), size(data,1));
yaxis = logspace(log10(lowlim), log10(hilim), size(data,2));

figure(h)
subplot(sampleCt,2,6)
surf(xaxis,yaxis,data./sum(sum(data)))
% colormap('gray')
% colormap(flipud(colormap))
shading flat
view([0 90])
set(gca,'XScale','log','YScale','log','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
line([lowlim hilim],[lowlim hilim],[max(max(data)) max(max(data))],'LineWidth',2,'Color','red','LineStyle','--')
xlabel({'$\it{T}\rm{_2 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
ylabel({'$\it{T}\rm{_1 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
title('Double CHIRP','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
colorbar
xlim(viewlim)
ylim(viewlim)

%% export

print -depsc2 '/Users/tyler/Dropbox/Manuscripts/CHIRP/ExcelTableData/AllT1T2Maps.eps'
