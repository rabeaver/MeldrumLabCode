clear
clc
close all
%%

data = load('C:\users\tkmeldrum\Dropbox\Data\CHIRP\Pretty2.out');

xaxis = logspace(-4, 0, size(data,1));
yaxis = logspace(-4, 0, size(data,2));

figure(2)
surf(xaxis,yaxis,data)
colormap('gray')
colormap(flipud(colormap))
shading flat
view([0 90])
set(gca,'XScale','log')
set(gca,'YScale','log')
line([1e-4 1],[1e-4 1],[max(max(data)) max(max(data))],'LineWidth',2,'Color','red','LineStyle','--')
xlabel('\it{T}\rm{_2 [s]}')
ylabel('\it{T}\rm{_1 [s]}')
