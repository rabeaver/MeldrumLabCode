clear
close all
clc


%%

dataDir = '/Users/tyler/Dropbox/Data/Perugia_May_2012/21May2012/IlPerugino/blue/2/';
parfilestem = 'acqu';

biexp_T2fit = 1;
T2rough_fit = 0;
omit_point = 1;
autoPhase = 1;
plotRaw = 1;


%%

%Load all data
[params,measDepth,depth,T1time,allT1Re,allT1Im,allT1Mag,T2time,allT2Re,allT2Im,allT1Phased,allT2Phased,T2rough,T2roughPhased,T2phaseAngle] = loadDataFromDepthInfo(dataDir,parfilestem,autoPhase);
% [params,measDepth,depth,T1time,allT1Re,allT1Im,allT1Mag,T2time,allT2Re,allT2Im,~,~,T2rough,~,~] = loadDataFromDepthInfo(dataDir,parfilestem,autoPhase);


% Plot all data
if plotRaw == 1;
%     figure(1)
%     hold all
%     h(1) = subplot(2,2,1);
%     meshc(measDepth,T1time,allT1Re)
%     xlabel('Depth (um)')
%     ylabel('time (ms)')
%     zlabel('signal intensity (arb)')
%     title('T1 Real')
%     h(2) = subplot(2,2,2);
%     meshc(measDepth,T1time,allT1Im)
%     xlabel('Depth (um)')
%     ylabel('time (ms)')
%     zlabel('signal intensity (arb)')
%     title('T1 Imag')
%     h(3) = subplot(2,2,3);
%     meshc(measDepth,T2time(omit_point+1:end),allT2Re(omit_point+1:end,:))
%     xlabel('Depth (um)')
%     ylabel('time (ms)')
%     zlabel('signal intensity (arb)')
%     title('T2 Real')
%     h(4) = subplot(2,2,4);
%     meshc(measDepth,T1time,allT1Mag)
%     xlabel('Depth (um)')
%     ylabel('time (ms)')
%     zlabel('signal intensity (arb)')
%     title('T1 Mag')
%     
    % Detailed T2 plot
    profileFig = figure('Color','none');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperType', 'A4');
    set(gcf, 'PaperOrientation', 'landscape');
    
    hold on
    plot(depth,real(T2roughPhased(omit_point+1,:)),'-k')
    line([min(depth);max(depth)],[0;0],'Color','k')
    
    for i=1:length(measDepth)
        xind = measDepth(i);
        yind = real(T2roughPhased(omit_point+1,find(depth == measDepth(i))));
        indlabel = num2str(i);
        
        yprev = real(T2roughPhased(omit_point+1,find(depth == measDepth(i))-1));
        yafter = real(T2roughPhased(omit_point+1,find(depth == measDepth(i))+1));
        
        if(yind < yprev && yind < yafter) %min
            text(xind,yind,indlabel,'VerticalAlignment','top','HorizontalAlignment','center');
        elseif(yind < yprev && yind > yafter) %pos slope
            text(xind,yind,indlabel,'VerticalAlignment','middle','HorizontalAlignment','right');
        elseif(yind > yprev && yind < yafter) %neg slope
            text(xind,yind,indlabel,'VerticalAlignment','middle','HorizontalAlignment','left');
        elseif(yind > yprev && yind > yafter) %max
            text(xind,yind,indlabel,'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
    end

    
%     figure(3)
%     mesh(depth,T2time(omit_point+1:end),real(T2roughPhased(omit_point+1:end,:)))
end


% Fit rough data to T2 decay
if T2rough_fit == 1
    guesses = [0;max(max(T2roughPhased));T2time(round(length(T2time)/8))];
    CI = 90; %desired confidence interval in percent
    ypred_T2_rough = zeros(1001,length(depth));
    
    for n = 1:1:length(depth)
        try
            [xfit_T2,ypred_T2_rough(:,n),T2_coeffs] = monodecay_t2fit(T2time(omit_point+1:end),real(T2roughPhased(omit_point+1:end,n)),guesses,CI);
            y0_T2_rough(n,:) = T2_coeffs(1,1:2);
            A_T2_rough(n,:) = T2_coeffs(2,1:2);
            tau_T2_rough(n,:) = T2_coeffs(3,1:2);
        catch
            y0_T2_rough(n,:) = [0;0];
            A_T2_rough(n,:) = [0;0];
            tau_T2_rough(n,:) = [0;0];
        end
    end
    tau_T2_rough(:,3) = tau_T2_rough(:,2)./tau_T2_rough(:,1)*100;
    
end

% Fit regular data to T2 decay
guesses = [0;max(max(allT2Re));T2time(round(length(T2time)/8))];
CI = 90; %desired confidence interval in percent
ypred_T2 = zeros(1001,size(measDepth,2));

for n = 1:1:size(measDepth,2)
    try
        [xfit_T2,ypred_T2(:,n),T2_coeffs,T2_resid(:,n)] = monodecay_t2fit(T2time(omit_point+1:end),real(allT2Phased(omit_point+1:end,n)),guesses,CI);
        y0_T2(n,:) = T2_coeffs(1,1:2);
        A_T2(n,:) = T2_coeffs(2,1:2);
        tau_T2(n,:) = T2_coeffs(3,1:2);
    catch
        y0_T2(n,:) = [0,0];
        A_T2(n,:) = [0,0];
        tau_T2(n,:) = [0,0];
        T2_resid(:,n) = zeros(length(T2time)-omit_point,1);
    end
end

tau_T2(:,3) = tau_T2(:,2)./tau_T2(:,1)*100;

% optional biexponential fit
if biexp_T2fit == 1;
    guesses = [0;max(max(allT2Re));T2time(round(length(T2time)/20));max(max(allT2Re))/4;T2time(round(length(T2time)/5))];
    CI = 90; %desired confidence interval in percent
    ypred_T2_bi = zeros(1001,size(measDepth,2));
    
    for n = 1:1:size(measDepth,2)
        try
            [xfit_T2,ypred_T2_bi(:,n),T2_coeffs_bi,T2_resid_bi(:,n)] = bidecay_t2fit(T2time(omit_point+1:end),real(allT2Phased(omit_point+1:end,n)),guesses,CI);
            y0_T2_bi(n,:) = T2_coeffs_bi(1,1:2);
            A_T2_bi(n,:) = T2_coeffs_bi(2,1:2);
            tau_T2_bi(n,:) = T2_coeffs_bi(3,1:2);
            A2_T2_bi(n,:) = T2_coeffs_bi(4,1:2);
            tau2_T2_bi(n,:) = T2_coeffs_bi(5,1:2);
        catch
            guesses = [y0_T2(n,1);A_T2(n,1);tau_T2(n,1);A_T2(n,1)/5;tau_T2(n,1)*5];
            try
                [xfit_T2,ypred_T2_bi(:,n),T2_coeffs_bi,T2_resid_bi(:,n)] = bidecay_t2fit(T2time(omit_point+1:end),real(allT2Phased(omit_point+1:end,n)),guesses,CI);
                y0_T2_bi(n,:) = T2_coeffs_bi(1,1:2);
                A_T2_bi(n,:) = T2_coeffs_bi(2,1:2);
                tau_T2_bi(n,:) = T2_coeffs_bi(3,1:2);
                A2_T2_bi(n,:) = T2_coeffs_bi(4,1:2);
                tau2_T2_bi(n,:) = T2_coeffs_bi(5,1:2);
            catch
                y0_T2_bi(n,:) = [0,0];
                A_T2_bi(n,:) = [0,0];
                tau_T2_bi(n,:) = [0,0];
                A2_T2_bi(n,:) = [0,0];
                tau2_T2_bi(n,:) = [0,0];
            end
        end
        
    end
    
    tau_T2_bi(:,3) = tau_T2_bi(:,2)./tau_T2_bi(:,1)*100;
    tau2_T2_bi(:,3) = tau2_T2_bi(:,2)./tau2_T2_bi(:,1)*100;
    
end



% Fit to T1 decay
guesses = [0;max(max(real(allT1Phased)));T1time(round(length(T1time)/6))];
CI = 90; %desired confidence interval in percent
ypred_T1 = zeros(1001,size(measDepth,2));

for n = 1:1:size(measDepth,2)
    [xfit_T1,ypred_T1(:,n),T1_coeffs,T1_resid(:,n)] = T1_fit(T1time,real(allT1Phased(:,n)),guesses,CI);
    y0_T1(n,:) = T1_coeffs(1,1:2);
    A_T1(n,:) = T1_coeffs(2,1:2);
    tau_T1(n,:) = T1_coeffs(3,1:2);
end

tau_T1(:,3) = tau_T1(:,2)./tau_T1(:,1)*100;

% Plot data and fits and select which to save

T1T2fig = figure('Units','normalized','OuterPosition',[0.15 0.15 0.7 0.7],'Color','none'); %,'MenuBar','none');
n = size(measDepth,2);
str{1} = '0';

for i = 1:n
    subplot(6,n,[i i+n])
    hold on
    %     plot(xfit_T1,ypred_T1(:,i))
    %     plot(T1time,allT1Re(:,i),'ob')
    %     plot(T1time,allT1Im(:,i),'or')
    %     plot(T1time,allT1Mag(:,i),'ob')
    plot(T1time,real(allT1Phased(:,i)),'ok')
    plot(xfit_T1,ypred_T1(:,i),'-k')
    ylim([-0.5*max(max(allT1Mag)) 1.1*max(max(allT1Mag))])
    xlim([0 max(T1time)])
    title(i)

    subplot(6,n,2*n+i)
    hold on
    plot(T1time,T1_resid(:,i),'-ok')
    line([0;max(T1time)],[0;0],'Color','k')
    ylim([-1.1*max(max(T1_resid)) 1.1*max(max(T1_resid))])
    xlim([0 max(T1time)])

    
    subplot(6,n,[3*n+i i+4*n])
    hold on
    %     plot(T2time,allT2Re(:,i),'-ob')
    %     plot(T2time,allT2Im(:,i),'-or')
    plot(T2time,real(allT2Phased(:,i)),'ok')
    plot(xfit_T2,ypred_T2(:,i),'-k')
    if biexp_T2fit == 1;
        plot(xfit_T2,ypred_T2_bi(:,i),'-.k')
    end
    %     plot(xfit_T2_lin,ypred_T2_lin(:,i),'-r')
    limits = max(max(real(allT2Phased)));
    ylim([-0.1*limits 1.1*limits])
    xlim([0 max(T2time)])

    subplot(6,n,5*n+i)
    hold on
    plot(T2time(omit_point + 1:end),T2_resid(:,i),'-k')
    plot(T2time(omit_point + 1:end),T2_resid_bi(:,i),'-.k')
    %     plot(T2time(omit_point + abs(linearizedFit-1):end),T2_resid(:,i),'-b')
    line([0;max(T2time)],[0;0],'Color','k')
    %     bar([T2time(omit_point+1:end),T2_resid(:,i),T2_resid_bi(:,i)])
    %     bar(T2time(omit_point+1:end),T2_resid_bi(:,i))
    ylim([-0.2*limits 0.2*limits])
    xlim([0 max(T2time)])
end