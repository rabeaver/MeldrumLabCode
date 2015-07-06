clear
clc
close all

%give file dir, file name (the *.out file from Prospa export2d), T2 and D limits (should be the same for Naproxed
%stuff), and number of points in inverted data.
datadir = '/Users/tyler/Desktop/';
datafile = '6_4_BSA_NPNa_T2D_Normalizes_Inverted.out';
T2lims = [1e-4 1e1];
Dlims = [1e-11 1e-9];
nPts = 750;

%this calculates the axes based on the limits above
Daxis = logspace(log10(Dlims(1)),log10(Dlims(2)),nPts);
T2axis = logspace(log10(T2lims(1)),log10(T2lims(2)),nPts);

%load the data and remove 0 values (replace with NaN)
data = load(strcat(datadir,datafile));
data(data == 0) = NaN; 

%this part requires teh extrema2.m and extrema.m functions, available from
%the mathworks file exchange. Finds the indices of the various peaks.
[xymax,smax,~,~] = extrema2(data);
[T2ind,Dind] = ind2sub([nPts,nPts],smax);

%plot the T2D data
figure(5)
surf(T2axis,Daxis,data)
shading flat
set(gca,'XScale','log','YScale','log')
xlim(T2lims)
ylim(Dlims)
ylabel('D [m^2 s^{-1}]')
xlabel('T_2 [s]')
view([0,90])

%for each peak present in the sample, make a contour line showing the 50%
%level. This countour is stored as "c", with some points designating the
%0.5 level and the actual height of that countour, and the rest defining
%the shape of the contour. The shape of the countour is plotted and the D_n
%and T2_n are arrays that hold the center value for the peak and the total 
%range for the 50% contour level as [low mean high].
%This is repeated for each peak present (up to three currently).
%Note: may need to adjust limits for drawing contours for peaks greater
%than 1--the figure will look weird if you don't. An explanation of how to
%do so is below.

%peak 1
n = 1;
figure(n)
[c,~] = contour(T2axis,Daxis,data./data(T2ind(n),Dind(n)),[0.5,0.5]);
plot(c(1,2:end),c(2,2:end),'LineWidth',2);
set(gca,'XScale','log','YScale','log')
xlim(T2lims)
ylim(Dlims)
ylabel('D [m^2 s^{-1}]')
xlabel('T_2 [s]')
D_1 =[ min(c(2,2:end)) Daxis(T2ind(n)) max(c(2,2:end))];
T2_1 = [ min(c(1,2:end)) T2axis(Dind(n)) max(c(1,2:end))];

%peak 2 if present
if length(T2ind) > 1
    n = 2;
    figure(n)
    [c,~] = contour(T2axis,Daxis,data./data(T2ind(n),Dind(n)),[0.5,0.5]);
    plot(c(1,2:382),c(2,2:382),'LineWidth',2);             % will need to update the number c(1,2:XXX) for different data sets. For secondary peaks, the countour line for the main peak still shows up in two places, so need to specifiy the end of the first peak. Just do this graphically.
    set(gca,'XScale','log','YScale','log')
    xlim(T2lims)
    ylim(Dlims)
    ylabel('D [m^2 s^{-1}]')
    xlabel('T_2 [s]')
    D_2 =[ min(c(2,2:382)) Daxis(T2ind(n)) max(c(2,2:382))];    % same here and below
    T2_2 = [ min(c(1,2:382)) T2axis(Dind(n)) max(c(1,2:382))];
end

%peak 3 if present
if length(T2ind) > 2
    n = 3;
    figure(n)
    [c,~] = contour(T2axis,Daxis,data./data(T2ind(n),Dind(n)),[0.5,0.5]);
    plot(c(1,2:382),c(2,2:382));
    set(gca,'XScale','log','YScale','log')
    xlim(T2lims)
    ylim(Dlims)
    ylabel('D [m^2 s^{-1}]')
    xlabel('T_2 [s]')
    D_2 =[ min(c(2,2:382)) Daxis(T2ind(n)) max(c(2,2:382))];
    T2_2 = [ min(c(1,2:382)) T2axis(Dind(n)) max(c(1,2:382))];
end


%%
% figure(1)
% subplot(2,2,1)
% pcolor(T2axis,Daxis,data)
% shading flat
% set(gca,'XScale','log','YScale','log')
% ylabel('D [m^2 s^{-1}]')
% xlabel('T_2 [s]')
% view([0,90])
% subplot(2,2,2)
% surf(T2axis,Daxis,data)
% shading flat
% set(gca,'XScale','log','YScale','log')
% ylabel('D [m^2 s^{-1}]')
% xlabel('T_2 [s]')
% view([-90,0])
% subplot(2,2,4)
% surf(T2axis,Daxis,data)
% shading flat
% set(gca,'XScale','log','YScale','log')
% ylabel('D [m^2 s^{-1}]')
% xlabel('T_2 [s]')
% view([0,0])


