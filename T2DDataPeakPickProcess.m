clear
clc
close all

%give file dir, file name (the *.out file from Prospa export2d), T2 and D limits (should be the same for Naproxed
%stuff), and number of points in inverted data.
datadir = '/Users/tyler/Desktop/BLANKS/';
datafile = 'T2_D_7_3_NPNa_Only_1_25JUN.out';
T2lims = [1e-4 1e0];
Dlims = [1e-11 1e-8];
contourLevel = 0.90;

%load the data and remove 0 values (replace with NaN)
data = load(strcat(datadir,datafile));
data(data == 0) = NaN; 
nPts = size(data,1);

%this calculates the axes based on the limits above
Daxis = logspace(log10(Dlims(1)),log10(Dlims(2)),nPts);
T2axis = logspace(log10(T2lims(1)),log10(T2lims(2)),nPts);


%this part requires the extrema2.m and extrema.m functions, available from
%the mathworks file exchange. Finds the indices of the various peaks.
[xymax,smax,~,~] = extrema2(data);
[T2ind,Dind] = ind2sub([nPts,nPts],smax);

% cm = colormap(gray);

% %plot the T2D data
figure(1)
surf(T2axis,Daxis,data)
colormap(flipud(gray));
shading flat
set(gca,'XScale','log','YScale','log','XTick', [1e-4; 1e-3; 1e-2; 1e-1; 1e0; 1e1],'FontSize',18)
xlim(T2lims)
ylim(Dlims)
ylabel('\itD\rm [m^2 s^{-1}]')
xlabel('\itT\rm_2 [s]')
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
for n = 1:length(T2ind); 
% n = 5;
    

    [c.n{n},~] = contour(T2axis,Daxis,data./data(T2ind(n),Dind(n)),[contourLevel,contourLevel]);
end

for n = 1:length(T2ind);
    figure(n)
    plot(c.n{n}(1,:))
end

%%
    ll = [2,   2,  90,  2, 296,  49,  79, 311]; %how to automate ll and mm?
    mm = [50, 25, 204, 88, 375, 181, 104, 596];
    lastPt = 2;

  
close all     
for n = 1:lastPt; 
%     figure(length(T2ind)+n+1)
%     plot(c.n{n}(1,ll(n):mm(n)),c.n{n}(2,ll(n):mm(n)),'LineWidth',2);             % will need to update the number c(1,2:XXX) for different data sets. For secondary peaks, the countour line for the main peak still shows up in two places, so need to specifiy the end of the first peak. Just do this graphically.
%     set(gca,'XScale','log','YScale','log')
%     xlim(T2lims)
%     ylim(Dlims)
%     ylabel('D [m^2 s^{-1}]')
%     xlabel('T_2 [s]'
    D(n,:) =[ min(c.n{n}(2,ll(n):mm(n))) Daxis(T2ind(n)) max(c.n{n}(2,ll(n):mm(n)))];
    T2(n,:) = [ min(c.n{n}(1,ll(n):mm(n))) T2axis(Dind(n)) max(c.n{n}(1,ll(n):mm(n)))]*1e3;
end

figure(1)
surf(T2axis,Daxis,data)
colormap(flipud(gray));
shading flat
set(gca,'XScale','log','YScale','log','XTick', [1e-4; 1e-3; 1e-2; 1e-1; 1e0; 1e1],'FontSize',18)
xlim(T2lims)
ylim(Dlims)
ylabel('\itD\rm [m^2 s^{-1}]')
xlabel('\itT\rm_2 [s]')
view([0,90])
hold on

for n = 1:lastPt;
    plot3(c.n{n}(1,ll(n):mm(n)),c.n{n}(2,ll(n):mm(n)),ones(1,mm(n)-ll(n)+1),'-r','LineWidth',3); 
    text(min(c.n{n}(1,ll(n):mm(n))), max(c.n{n}(2,ll(n):mm(n))),1,int2str(n));
end
set(gca,'XScale','log','YScale','log')
    xlim(T2lims)
    ylim(Dlims)
    ylabel('D [m^2 s^{-1}]')
    xlabel('T_2 [s]')

%plot the T2D data
% figure
% surf(T2axis,Daxis,data)
% colormap(flipud(gray));
% shading flat
% set(gca,'XScale','log','YScale','log','XTick', [1e-4; 1e-3; 1e-2; 1e-1; 1e0; 1e1],'FontSize',18)
% xlim(T2lims)
% ylim(Dlims)
% ylabel('\itD\rm [m^2 s^{-1}]')
% xlabel('\itT\rm_2 [s]')
% view([0,90])

%%
%peak 2 if present
if length(T2ind) > 1
    n = 2;
    ll2 = 1;
    mm2 = 5;
    figure(n)
    [c2,~] = contour(T2axis,Daxis,data./data(T2ind(n),Dind(n)),[contourLevel,contourLevel]);
    plot(c2(1,ll2:mm2),c2(2,ll2:mm2),'LineWidth',2);             % will need to update the number c(1,2:XXX) for different data sets. For secondary peaks, the countour line for the main peak still shows up in two places, so need to specifiy the end of the first peak. Just do this graphically.
    set(gca,'XScale','log','YScale','log')
    xlim(T2lims)
    ylim(Dlims)
    ylabel('D [m^2 s^{-1}]')
    xlabel('T_2 [s]')
    D(n,:) =[ min(c2(2,ll2:mm2)) Daxis(T2ind(n)) max(c2(2,ll2:mm2))];    % same here and below
    T2(n,:) = [ min(c2(1,ll2:mm2)) T2axis(Dind(n)) max(c2(1,ll2:mm2))]*1e3;
end


figure(5)
hold on
plot3(c1(1,ll1:mm1),c1(2,ll1:mm1),ones(1,mm1-ll1+1),'-r','LineWidth',3); 
plot3(c2(1,ll2:mm2),c2(2,ll2:mm2),ones(1,mm2-ll2+1),'-r','LineWidth',3); 

%%
%peak 3 if present
if length(T2ind) > 2
    n = 3;
    ll3 = 17;
    mm3 = 47;
    figure(n)
    [c3,~] = contour(T2axis,Daxis,data./data(T2ind(n),Dind(n)),[contourLevel,contourLevel]);
    plot(c3(1,ll3:mm3),c3(2,ll3:mm3),'LineWidth',2);             % will need to update the number c(1,2:XXX) for different data sets. For secondary peaks, the countour line for the main peak still shows up in two places, so need to specifiy the end of the first peak. Just do this graphically.
    set(gca,'XScale','log','YScale','log')
    xlim(T2lims)
    ylim(Dlims)
    ylabel('D [m^2 s^{-1}]')
    xlabel('T_2 [s]')
    D(n,:) =[ min(c3(2,ll3:mm3)) Daxis(T2ind(n)) max(c3(2,ll3:mm3))];    % same here and below
    T2(n,:) = [ min(c3(1,ll3:mm3)) T2axis(Dind(n)) max(c3(1,ll3:mm3))]*1e3;
end

figure(5)
hold on
plot3(c1(1,ll1:mm1),c1(2,ll1:mm1),ones(1,mm1-ll1+1),'-r','LineWidth',3); 
plot3(c2(1,ll2:mm2),c2(2,ll2:mm2),ones(1,mm2-ll2+1),'-r','LineWidth',3); 
plot3(c3(1,ll3:mm3),c3(2,ll3:mm3),ones(1,mm3-ll3+1),'-r','LineWidth',3); 

%%
%peak 4 if present
if length(T2ind) > 3
    n = 4;
    ll4 = 160;
    mm4 = 212;
    figure(n)
    [c4,~] = contour(T2axis,Daxis,data./data(T2ind(n),Dind(n)),[contourLevel,contourLevel]);
    plot(c4(1,ll4:mm4),c4(2,ll4:mm4),'LineWidth',2);             % will need to update the number c(1,2:XXX) for different data sets. For secondary peaks, the countour line for the main peak still shows up in two places, so need to specifiy the end of the first peak. Just do this graphically.
    set(gca,'XScale','log','YScale','log')
    xlim(T2lims)
    ylim(Dlims)
    ylabel('D [m^2 s^{-1}]')
    xlabel('T_2 [s]')
    D(n,:) =[ min(c4(2,ll4:mm4)) Daxis(T2ind(n)) max(c4(2,ll4:mm4))];    % same here and below
    T2(n,:) = [ min(c4(1,ll4:mm4)) T2axis(Dind(n)) max(c4(1,ll4:mm4))]*1e3;
end

figure(5)
hold on
plot3(c1(1,ll1:mm1),c1(2,ll1:mm1),ones(1,mm1-ll1+1),'-r','LineWidth',3); 
plot3(c2(1,ll2:mm2),c2(2,ll2:mm2),ones(1,mm2-ll2+1),'-r','LineWidth',3); 
plot3(c3(1,ll3:mm3),c3(2,ll3:mm3),ones(1,mm3-ll3+1),'-r','LineWidth',3); 
plot3(c4(1,ll4:mm4),c4(2,ll4:mm4),ones(1,mm4-ll4+1),'-r','LineWidth',3); 


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


