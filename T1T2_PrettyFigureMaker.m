clear
clc
close all

filedir = '/Users/tyler/Dropbox/Data/CHIRP/Dec2015/ConcentricSamples_01Dec2015/';
filename = 'Double_CHIRP_1024_250PCD_08Dec2015_result';
data = load(strcat(filedir,'out files/ ',filename,'.out'));
lowlim = 1e-4;
hilim = 1e-1;

viewLim = [lowlim hilim];

xaxis = logspace(log10(lowlim), log10(hilim), size(data,1));
yaxis = logspace(log10(lowlim), log10(hilim), size(data,2));

figure
pcolor(xaxis,yaxis,data./sum(sum(data)))
set(gca,'XScale','log','YScale','log','TickLabelInterpreter', 'latex','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
set(gca,'XTickLabel','','YTickLabel','');
colormap('gray')
colormap(flipud(colormap))
shading flat
% line(viewLim,viewLim,'LineWidth',2,'Color','red','LineStyle','--')
% xlabel({'$\it{T}\rm{_2 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
% ylabel({'$\it{T}\rm{_1 [s]}$'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times')
% title('15 mM Gd^{3+}(aq), Traditional','FontUnits','points','FontWeight','normal','FontSize',9,'FontName','Times')
% colorbar
xlim(viewLim)
ylim(viewLim)
set(gca,'LooseInset',get(gca,'TightInset'))

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize', [5 5]);
set(gcf,'PaperPosition',[0 0 5 5]);
set(gcf,'PaperPositionMode','Manual');

% print(gcf, '-depsc2',  strcat(filedir, 'figures/', filename,'.eps'));
print(gcf, '-dtiff', '-r300', strcat(filedir, 'figures/', filename,'.tif'));
