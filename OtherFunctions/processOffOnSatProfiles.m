function [spectrum, FID, ph0, contrast, onresfit, tauon, offresfit, tauoff, contrastfit, taucontrast] = processOffOnSatProfiles(lsfid,lb,baselinepts,datapath,fidfolder)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user defined parameters for processing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bubble = readpar(strcat(datapath,fidfolder,'/procpar'),'bubble '); %bubb;e time
settle = readpar(strcat(datapath,fidfolder,'/procpar'),'settle '); %settle time
p1 = readpar(strcat(datapath,fidfolder,'/procpar'),'p1 '); %pulse time
satcycle = readpar(strcat(datapath,fidfolder,'/procpar'),'satcycle '); %satcycles
satpwr = readpar(strcat(datapath,fidfolder,'/procpar'),'satpwr '); %saturation power
satfreq = readpar(strcat(datapath,fidfolder,'/procpar'),'satfreq '); %saturation frequency
at = readpar(strcat(datapath,fidfolder,'/procpar'),'at '); %acquisition time
sw = readpar(strcat(datapath,fidfolder,'/procpar'),'sw '); %spectral width
nt = readpar(strcat(datapath,fidfolder,'/procpar'),'nt '); %number of averages

sattime = p1 * satcycle;

np = 18; %max(size(satfreq))*max(size(sattime));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process FID into spectrum, FID, and ph0 information %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[spectrum, FID, ph0, rawFID] = processFID(fidfolder,datapath,lsfid,baselinepts,lb,1,sw,at,1,6000,8000);

integrationRange = 500;
[C,I] = max(abs(spectrum(:,:)));
for i=2:2:np
       if I(i-1) <= integrationRange
           I(i-1) = integrationRange + 1;
       end
       if I(i) - integrationRange < 1
           I(i) = integrationRange + 1;
       end
    data4int(i-1) = squeeze(sum(abs(spectrum(I(i-1) - integrationRange:I(i-1) + integrationRange,i-1))));
    data4int(i) = squeeze(sum(abs(spectrum(I(i-1) - integrationRange:I(i-1) + integrationRange,i))));
end
%data4=squeeze(sum(real(spectrum((I(1)-200):(I(1) + 200),:)))); %integrate the peak area +/- 250 pts from the highest point
data4 = data4int';

xaxis = 1:1:size(spectrum(:,1));
xaxis = xaxis';


% figure
% for i = 1:1:21
%     subplot(5,5,i)
%     plot(real(spectrum(:,i)))
% end

%% Plot of Contrast
    
    for i=2:2:np
        contrast(i/2) = 1- (data4(i-1) - data4(i))/(data4(i-1));
        time(i/2) = sattime(i/2);
        offres(i/2) = data4(i - 1);
        onres(i/2) = data4(i);
    end
    
%     figure
%     hold on
%     plot(time,offres,'-k')
%     plot(time,onres,'-r')
%     
%     figure
%     plot(time,contrast);
%     
    
%% Fitting to Exponential Decay
disp('On Res Saturation')
[onresfit,tauon] = fittingRoutineNoSave(onres,time,.90);
disp('Off Res Saturation')
[offresfit,tauoff] = fittingRoutineNoSave(offres,time,.90);
disp('Contrast')
[contrastfit,taucontrast] = fittingRoutineTotalExpNoSave(contrast,time,.90);

% Final Plotting
figure
    subplot(2,1,1)
%         title('Saturation Profile for 12dB cw')
        hold on
        plot(time,onresfit,'-b')
        plot(time,offresfit,'-k')
        ylabel('normalized signal (arb)')
        legend('On res','Off res')
        ylim([0 1])
        hold off

    subplot(2,1,2)
        hold on
        plot(time,-contrastfit + 1,'-r')
        plot(time,-contrast + 1,':r')
        hold off
        ylim([0 1])
        xlabel('saturation time (s)')
        ylabel('contrast [ (off-on)/off) ]')
        legend('Fit','Experimental','Location','Best')