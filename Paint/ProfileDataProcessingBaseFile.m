% This is my new base file for data processing.

clear
clc
close all
% try
%     matlabpool
% catch
% end

%%

info(1).parfilename = 'acqu';
info(2).parfilename = 'acqu';
info(1).datafilename = 'VarnishThick_10_3';
info(2).datafilename = 'VarnishThick_10_3';
info(1).dirstem = '/Users/jaredking/Documents/Research Files and Data/Paint/Varnish Model Data/Dry Varnish Scans/ThickVarnish/VarnishThick_10_3/';
info(2).dirstem = '/Users/jaredking/Documents/Research Files and Data/Paint/Varnish Model Data/Dry Varnish Scans/ThickVarnish/VarnishThick_10_3/';
omitpoints = 0;

fitopts = statset('MaxIter',5000,'TolX',1e-14,'UseParallel',true,'Display','off');
% guesses = [0; 0.01; 100];% 0.1; 2];
guesses = [0.1; 100];% 0.1; 2];
guessesbi = [0.01; 500; 0.15; 4000];
% lb = [0;0;0];%0;0];
% ub = [inf,inf,inf];%,inf,inf];


for k = 1:2
dir = strcat(info(k).dirstem,num2str(1),'/');
cd(dir);

params.acqTime(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'acqTime');
params.bandwidth(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'bandwidth');
params.nrScans(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'nrScans');
params.rxPhase(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'rxPhase');
params.rxGain(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'rxGain');
params.nrPts(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'nrPnts');
params.repTime(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'repTime');
params.b1Freq(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'b1Freq');
params.nrEchoes(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'nrEchoes');
params.echoTime(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'echoTime');
params.nrExp(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'nrExp');
params.initDepth(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'initDepth');
params.finalDepth(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'finalDepth');
params.stepSize(k) = readpar_Kea(strcat(info(k).parfilename,'.par'),'stepSize');

% if params.nrExp > 1;
%     parfor i= 1:params.nrExp
%         numID = num2str(i);
%         real_data = load(strcat(filename,'0',numID,'-decaysRe.dat'));
%         imag_data = load(strcat(filename,'0',numID,'-decaysIm.dat'));
%         temp_data = complex(real_data(omitpoints+1:end,2:end),imag_data(omitpoints+1:end,2:end));
%         data(:,:,i) = real(autophase(temp_data,1));
%         contrast_data = load(strcat(filename,'1.dat'));
%     end
% else


    real_data = load(strcat(info(k).datafilename,'-decaysRe.dat'));
    imag_data = load(strcat(info(k).datafilename,'-decaysIm.dat'));
    temp_data = complex(real_data(omitpoints+1:end,2:end),imag_data(omitpoints+1:end,2:end));
    data = real(autophase(temp_data,1));
    contrast_data = load(strcat(info(k).datafilename,'.dat'));
    
% end

    final_data(k).position = contrast_data(:,1);
    final_data(k).echoVector = params.echoTime(k)*(omitpoints+1):params.echoTime(k):params.echoTime(k)*params.nrEchoes(k);
    final_data(k).data = data;
    
    xfit = 0:1:max(final_data(1).echoVector);
    
%     if k ~= 4
%         final_data(k).position = final_data(k).position(3:43);
%         final_data(k).data = final_data(k).data(:,3:43);
% %     elseif k = 2 || k = 3
% %         final_data(k).position = final_data(k).position(3:43);
% %         final_data(k).data = final_data(k).data(:,3:43);
%     end
    
   if k == 1; 
       for i=1:length(final_data(k).position)
           try
               [final_data(k).beta(:,i),final_data(k).residual(:,i),final_data(k).J(:,:,i)] = nlinfit(final_data(k).echoVector',final_data(k).data(:,i),@t2monofit_simple,guesses,fitopts);
%                final_data(k).beta(3,i) = 0;
%                final_data(k).beta(4,i) = 0;
           catch pm
               final_data(k).beta(:,i) = [0;0];%;0;0];
               final_data(k).J(:,:,i) = zeros(size(final_data(k).data,1),2);
               final_data(k).residual(:,i) = zeros(size(final_data(k).data,1),1);
           end
           [final_data(k).ypred(:,i),final_data(k).delta(:,i)] = nlpredci(@t2monofit_simple,xfit,final_data(k).beta(:,i),final_data(k).residual(:,i),'Jacobian',final_data(k).J(:,:,i));
           tt = nlparci(final_data(k).beta(:,i),final_data(k).residual(:,i),'jacobian',final_data(k).J(:,:,i));
           final_data(k).ci(:,i) = tt(:,2)-final_data(k).beta(:,i);
           [final_data(k).yvals(:,i),final_data(k).yvalsdelta(:,i)] = nlpredci(@t2monofit_simple,final_data(k).echoVector,final_data(k).beta(:,i),final_data(k).residual(:,i),'Jacobian',final_data(k).J(:,:,i));
       end
   elseif k ==2;
       for i=1:length(final_data(k).position)
        try
            [final_data(k).beta(:,i),final_data(k).residual(:,i),final_data(k).J(:,:,i)] = nlinfit(final_data(k).echoVector',final_data(k).data(:,i),@t2bifit_simple,guessesbi,fitopts);
            
        catch pm
            final_data(k).beta(:,i) = [0;0;0;0];
            final_data(k).J(:,:,i) = zeros(size(final_data(k).data,1),4);
            final_data(k).residual(:,i) = zeros(size(final_data(k).data,1),1);
        end
        [final_data(k).ypred(:,i),final_data(k).delta(:,i)] = nlpredci(@t2bifit_simple,xfit,final_data(k).beta(:,i),final_data(k).residual(:,i),'Jacobian',final_data(k).J(:,:,i));
        tt = nlparci(final_data(k).beta(:,i),final_data(k).residual(:,i),'jacobian',final_data(k).J(:,:,i));
        final_data(k).ci(:,i) = tt(:,2)-final_data(k).beta(:,i);
        [final_data(k).yvals(:,i),final_data(k).yvalsdelta(:,i)] = nlpredci(@t2bifit_simple,final_data(k).echoVector,final_data(k).beta(:,i),final_data(k).residual(:,i),'Jacobian',final_data(k).J(:,:,i));
       end
%    elseif k==2
%     guesses=[0.15; 10000; 0.15; 4000; 0.15];
%     lbs = [0;0;0;0;0];
%     ubs = [0.25;inf;0.25;inf;0.5];
% 
%     options = optimoptions('lsqcurvefit','Jacobian','on');
% 
%     for i = 1:length(final_data(k).position)
% 
%         try
%             [final_data(k).beta(:,i),final_data(k).residual(:,i),~,~,~,final_data(k).J(:,:,i)] = lsqcurvefit(@t2trifit,guesses,final_data(k).echoVector',final_data(k).data(:,i),lbs,ubs);
%         catch pm
%             pm
%             final_data(k).beta(:,i) = 0;
%             final_data(k).J(:,:,i) = zeros(size(final_data(k).data,1),6);
%             final_data(k).residual(:,i) = zeros(size(final_data(k).data,1),1);
%         end
% %         [final_data(k).ypred(:,i),final_data(k).delta(:,i)] = nlpredci(@t2trifit,xfit,final_data(k).beta(:,i),final_data(k).residual(:,i),'Jacobian',final_data(k).J(:,:,i));
% %         tt = nlparci(final_data(k).beta(:,i),final_data(k).residual(:,i),'jacobian',final_data(k).J(:,:,i));
% %         final_data(k).ci(:,i) = tt(:,2)-final_data(k).beta(:,i);
% %         [final_data(k).yvals(:,i),final_data(k).yvalsdelta(:,i)] = nlpredci(@t2trifit_simple,final_data(k).echoVector,final_data(k).beta(:,i),final_data(k).residual(:,i),'Jacobian',final_data(k).J(:,:,i));
%     end
   end
end   
%         final_data(k).yvalsrelerror = final_data(k).yvalsdelta./final_data(k).yvals;
% position(:,k) = contrast_data(:,1);

%     echoVector(:,k) = params.echoTime(k)*(omitpoints+1):params.echoTime(k):params.echoTime(k)*params.nrEchoes(k);
%     final_data = [echoVector(:,k),data];

%%
m = 2; %data to subtract from
l = 1; %data that is subtracted
guesses = [0.15; 3000];
guessesbi = [0.1; 100; 0.01; 1000];
for i = 1:(length(final_data(m).position)-2)
    subdata.data(:,i) = final_data(m).data(:,i+2) - final_data(l).yvals(:,i);
    try
        [subdata.beta(:,i),subdata.residual(:,i),subdata.J(:,:,i)] = nlinfit(final_data(m).echoVector',subdata.data(:,i),@t2monofit_simple,guesses,fitopts);
            
        catch pm
            subdata.beta(:,i) = [0;0];%;0;0];
            subdata.J(:,:,i) = zeros(size(final_data(m).data,1),2);
            subdata.residual(:,i) = zeros(size(final_data(m).data,1),1);
    end
        [subdata.ypred(:,i),subdata.delta(:,i)] = nlpredci(@t2monofit_simple,xfit,subdata.beta(:,i),subdata.residual(:,i),'Jacobian',subdata.J(:,:,i));
        tt = nlparci(subdata.beta(:,i),subdata.residual(:,i),'jacobian',subdata.J(:,:,i));
        subdata.ci(:,i) = tt(:,2)-subdata.beta(:,i);
end

%% So this is some weird new curve fitting tool that Tyler told me to use for a sample with weird stuff
m = 2; %data to subtract from
l = 1; %data that is subtracted
guesses=[0.15; 1e4; 0.15; 4000];
lbs = [0;4000;0;0];
ubs = [.25;inf;.25;4000];

options = optimoptions('lsqcurvefit','Jacobian','on');

for i = 1:length(final_data(m).position)
    subdata.data(:,i) = final_data(m).data(:,i) - final_data(l).yvals(:,i);
    try
        [subdata.beta(:,i),subdata.residual(:,i),~,~,~,subdata.J(:,:,i)] = lsqcurvefit(@t2bifit_simple,guesses,final_data(m).echoVector',subdata.data(:,i),lbs,ubs);
    catch
        subdata.beta(:,i) = [0;0;0;0];
        subdata.J(:,:,i) = zeros(size(final_data(m).data,1),4);
        subdata.residual(:,i) = zeros(size(final_data(m).data,1),1);
    end
end


%%

% for i = 1:length(final_data(m).position)
%     figure(i)
%     hold on
%     plot(final_data(m).echoVector,final_data(m).data(:,i),'xb')
%     plot(xfit,final_data(m).ypred(:,i),'--b')
%     plot(final_data(l).echoVector,final_data(l).data(:,i),'xk')
%     plot(xfit,final_data(l).ypred(:,i),'-k')
%     plot(final_data(m).echoVector,subdata.data(:,i),'xr')
%     plot(final_data(m).echoVector,subdata.data(:,i)-final_data(m).yvalsdelta(:,i),'--r')
%     plot(final_data(m).echoVector,subdata.data(:,i)+final_data(m).yvalsdelta(:,i),'--r')
%     plot(xfit,subdata.ypred(:,i),'-g')
%     plot(xfit,subdata.ypred(:,i)+subdata.delta(:,i),'--g')
%     plot(xfit,subdata.ypred(:,i)-subdata.delta(:,i),'--g')
%     title(strcat(num2str(final_data(m).position(i)),' um, No. ',num2str(i)))
%     xlim([0 max(final_data(m).echoVector)])
%     ylim([-0.01 0.2])
%     line([0 max(final_data(m).echoVector)],[0 0],'Color',[0 0 0])
% end

figure(1)
for i = 1:length(final_data(m).position)
    subplot(9,8,i)
    hold on
    plot(final_data(m).echoVector,final_data(m).data(:,i),'xb')
    plot(xfit,final_data(m).ypred(:,i),'--b')
    plot(final_data(l).echoVector,final_data(l).data(:,i),'xk')
    plot(xfit,final_data(l).ypred(:,i),'-k')
    plot(final_data(m).echoVector,subdata.data(:,i),'xr')
    plot(final_data(m).echoVector,subdata.data(:,i)-final_data(m).yvalsdelta(:,i),'--r')
    plot(final_data(m).echoVector,subdata.data(:,i)+final_data(m).yvalsdelta(:,i),'--r')
    plot(xfit,subdata.ypred(:,i),'-g')
    plot(xfit,subdata.ypred(:,i)+subdata.delta(:,i),'--g')
    plot(xfit,subdata.ypred(:,i)-subdata.delta(:,i),'--g')
    title(strcat(num2str(final_data(m).position(i)),' um'))
    xlim([0 max(final_data(m).echoVector)])
    ylim([-0.01 0.2])
    line([0 max(final_data(m).echoVector)],[0 0],'Color',[0 0 0])
end

%%
figure(2)
for i = 1:length(final_data(m).position)
    subplot(8,6,i)
    hold on
    plot(final_data(m).echoVector,subdata.residual(:,i),'-b')
    plot(final_data(m).echoVector,final_data(l).residual(:,i)+0.02,'-r')
    plot(final_data(m).echoVector,final_data(m).residual(:,i)-0.02,'-k')
%     plot(xfit,final_data(m).ypred(:,i),'--b')
%     plot(final_data(l).echoVector,final_data(l).data(:,i),'.k')
%     plot(xfit,final_data(l).ypred(:,i),'-k')
%     plot(final_data(m).echoVector,subdata.data(:,i),'.r')
%     plot(final_data(m).echoVector,subdata.data(:,i)-final_data(m).yvalsdelta(:,i),'--r')
%     plot(final_data(m).echoVector,subdata.data(:,i)+final_data(m).yvalsdelta(:,i),'--r')
%     plot(xfit,subdata.ypred(:,i),'-g')
%     plot(xfit,subdata.ypred(:,i)+subdata.delta(:,i),'--g')
%     plot(xfit,subdata.ypred(:,i)-subdata.delta(:,i),'--g')
    title(strcat(num2str(final_data(m).position(i)),' um'))
    xlim([0 max(final_data(m).echoVector)])
    ylim([-0.04 0.04])
    line([0 max(final_data(m).echoVector)],[0 0],'Color',[0 0 0])
    if i == length(final_data(m).position)
    legend('subtraction','dry','wet')
    end
end

%%
positionlims = [1 44];
figure(3)
subplot(2,3,1)
hold on
plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(2,positionlims(1):positionlims(2)),'-r')
plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(2,positionlims(1):positionlims(2))+final_data(1).ci(2,positionlims(1):positionlims(2)),'--r')
plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(2,positionlims(1):positionlims(2))-final_data(1).ci(2,positionlims(1):positionlims(2)),'--r')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(2,positionlims(1):positionlims(2)),'-k')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(2,positionlims(1):positionlims(2))+final_data(2).ci(2,positionlims(1):positionlims(2)),'--k')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(2,positionlims(1):positionlims(2))-final_data(2).ci(2,positionlims(1):positionlims(2)),'--k')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(4,positionlims(1):positionlims(2)),'-g')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(4,positionlims(1):positionlims(2))+final_data(2).ci(4,positionlims(1):positionlims(2)),'--g')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(4,positionlims(1):positionlims(2))-final_data(2).ci(4,positionlims(1):positionlims(2)),'--g')
% plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(2,positionlims(1):positionlims(2)))
% plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(2,positionlims(1):positionlims(2))+subdata.ci(2,positionlims(1):positionlims(2)),'--b')
% plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(2,positionlims(1):positionlims(2))-subdata.ci(2,positionlims(1):positionlims(2)),'--b')
ylim([0 3500])
title('T_2 vs. position')
ylabel('T_2 (us) from fit to subtracted data')
subplot(2,3,2)
hold on
% plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(2,positionlims(1):positionlims(2)),'-r')
% plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(2,positionlims(1):positionlims(2))+final_data(1).ci(2,positionlims(1):positionlims(2)),'--r')
% plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(2,positionlims(1):positionlims(2))-final_data(1).ci(2,positionlims(1):positionlims(2)),'--r')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(2,positionlims(1):positionlims(2)),'-k')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(2,positionlims(1):positionlims(2))+final_data(2).ci(2,positionlims(1):positionlims(2)),'--k')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(2,positionlims(1):positionlims(2))-final_data(2).ci(2,positionlims(1):positionlims(2)),'--k')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(4,positionlims(1):positionlims(2)),'-g')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(4,positionlims(1):positionlims(2))+final_data(2).ci(4,positionlims(1):positionlims(2)),'--g')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(4,positionlims(1):positionlims(2))-final_data(2).ci(4,positionlims(1):positionlims(2)),'--g')
% plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(2,positionlims(1):positionlims(2)))
% plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(2,positionlims(1):positionlims(2))+subdata.ci(2,positionlims(1):positionlims(2)),'--b')
% plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(2,positionlims(1):positionlims(2))-subdata.ci(2,positionlims(1):positionlims(2)),'--b')
ylim([0 20000])
title('T_2 vs. position')
ylabel('T_2 (us) from fit to subtracted data')
subplot(2,3,3)
hold on
% plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(2,positionlims(1):positionlims(2)),'-r')
% plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(2,positionlims(1):positionlims(2))+final_data(1).ci(2,positionlims(1):positionlims(2)),'--r')
% plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(2,positionlims(1):positionlims(2))-final_data(1).ci(2,positionlims(1):positionlims(2)),'--r')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(2,positionlims(1):positionlims(2)),'-k')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(2,positionlims(1):positionlims(2))+final_data(2).ci(2,positionlims(1):positionlims(2)),'--k')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(2,positionlims(1):positionlims(2))-final_data(2).ci(2,positionlims(1):positionlims(2)),'--k')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(4,positionlims(1):positionlims(2)),'-g')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(4,positionlims(1):positionlims(2))+final_data(2).ci(4,positionlims(1):positionlims(2)),'--g')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(4,positionlims(1):positionlims(2))-final_data(2).ci(4,positionlims(1):positionlims(2)),'--g')
plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(2,positionlims(1):positionlims(2)))
% plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(4,positionlims(1):positionlims(2)),'r')
plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(2,positionlims(1):positionlims(2))+subdata.ci(2,positionlims(1):positionlims(2)),'--b')
plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(2,positionlims(1):positionlims(2))-subdata.ci(2,positionlims(1):positionlims(2)),'--b')
ylim([0 15000])
title('T_2 vs. position')
ylabel('T_2 (us) from fit to subtracted data')
subplot(2,3,4)
hold on
plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(1,positionlims(1):positionlims(2)),'-r')
plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(1,positionlims(1):positionlims(2))+final_data(1).ci(1,positionlims(1):positionlims(2)),'--r')
plot(final_data(1).position(positionlims(1):positionlims(2)),final_data(1).beta(1,positionlims(1):positionlims(2))-final_data(1).ci(1,positionlims(1):positionlims(2)),'--r')
ylim([-0.05 0.3])
title('Amplitude vs. position, Dry')
xlabel('canvas <--- position (um) ---> paint')
ylabel('Amplitude (arb) from fit to subtracted data')
subplot(2,3,5)
hold on
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(1,positionlims(1):positionlims(2)),'-g')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(1,positionlims(1):positionlims(2))+final_data(2).ci(1,positionlims(1):positionlims(2)),'--g')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(1,positionlims(1):positionlims(2))-final_data(2).ci(1,positionlims(1):positionlims(2)),'--g')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(3,positionlims(1):positionlims(2)),'-k')
% plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(5,positionlims(1):positionlims(2)),'-r')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(3,positionlims(1):positionlims(2))+final_data(2).ci(3,positionlims(1):positionlims(2)),'--k')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(3,positionlims(1):positionlims(2))-final_data(2).ci(3,positionlims(1):positionlims(2)),'--k')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(1,positionlims(1):positionlims(2))+final_data(2).beta(3,positionlims(1):positionlims(2)),'-y')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(1,positionlims(1):positionlims(2))+final_data(2).beta(3,positionlims(1):positionlims(2))+(final_data(2).ci(1,positionlims(1):positionlims(2)).^2+final_data(2).ci(3,positionlims(1):positionlims(2)).^2).^0.5,'--y')
plot(final_data(2).position(positionlims(1):positionlims(2)),final_data(2).beta(1,positionlims(1):positionlims(2))+final_data(2).beta(3,positionlims(1):positionlims(2))-(final_data(2).ci(1,positionlims(1):positionlims(2)).^2+final_data(2).ci(3,positionlims(1):positionlims(2)).^2).^0.5,'--y')
ylim([-0.05 0.3])
title('Amplitude vs. position, Wet')
xlabel('canvas <--- position (um) ---> paint')
ylabel('Amplitude (arb) from fit to subtracted data')
subplot(2,3,6)
hold on
plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(1,positionlims(1):positionlims(2)))
% plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(3,positionlims(1):positionlims(2)),'r')
plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(1,positionlims(1):positionlims(2))+subdata.ci(1,positionlims(1):positionlims(2)),'--b')
plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(1,positionlims(1):positionlims(2))-subdata.ci(1,positionlims(1):positionlims(2)),'--b')
% plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(3,positionlims(1):positionlims(2)),'-m')
% plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(3,positionlims(1):positionlims(2))+subdata.ci(3,positionlims(1):positionlims(2)),'--m')
% plot(final_data(1).position(positionlims(1):positionlims(2)),subdata.beta(3,positionlims(1):positionlims(2))-subdata.ci(3,positionlims(1):positionlims(2)),'--m')
% plot(final_data(2).position(positionlims(1):positionlims(2)),subdata.beta(1,positionlims(1):positionlims(2))+subdata.beta(3,positionlims(1):positionlims(2)),'-c')
% plot(final_data(2).position(positionlims(1):positionlims(2)),subdata.beta(1,positionlims(1):positionlims(2))+subdata.beta(3,positionlims(1):positionlims(2))+(subdata.ci(1,positionlims(1):positionlims(2)).^2+subdata.ci(3,positionlims(1):positionlims(2)).^2).^0.5,'--c')
% plot(final_data(2).position(positionlims(1):positionlims(2)),subdata.beta(1,positionlims(1):positionlims(2))+subdata.beta(3,positionlims(1):positionlims(2))-(subdata.ci(1,positionlims(1):positionlims(2)).^2+subdata.ci(3,positionlims(1):positionlims(2)).^2).^0.5,'--c')
ylim([-0.05 0.3])
title('Amplitude vs. position, SubData')
xlabel('canvas <--- position (um) ---> paint')
ylabel('Amplitude (arb) from fit to subtracted data')

% figure(4)
% surf(final_data(1).position(positionlims(1):positionlims(2)),xfit,subdata.ypred(:,positionlims(1):positionlims(2)))
% shading interp
% zlim([0 0.2])
% caxis([0 0.2])

%%

figure(5)
for i = 1:length(final_data(m).position)
    subplot(9,8,i)
    hold on
    plot(final_data(1).echoVector,final_data(1).data(:,i),'-r')
    plot(xfit,final_data(1).ypred(:,i),'--r')
    plot(final_data(2).echoVector,subdata.data(:,i),'-k')
    plot(xfit,subdata.ypred(:,i),'--k')
    ylim([-0.01 0.1])
%     xlim([0 1000])
end
%%
% figure(1)
% hold on
% for i = 1:length(position)
%     subplot(6,7,i)
%     hold on
%     plot(xfit,t2monofit(beta(:,i),xfit),'-r')
%     plot(xfit,ypred(:,i),'-k')
%     plot(xfit,ypred(:,i)+delta(:,i),':b')
%     plot(xfit,ypred(:,i)-delta(:,i),':b')    
%     plot(sorted_data(:,1),sorted_data(:,i+1),'+b')
%     title(strcat(num2str(position(i)),' um'))
%     xlim([0 max(echoVector)])
%     ylim([-0.01 0.15])
% end
% 
% 
% figure(3)
% hold on
% plot(position,(beta(2,:)))
% ylim([0 0.25])
% xlim([1150 1650])
% 
% % figure(5)
% % hold on
% % plot(position,resnorm)
% 
% figure(5)
% hold on
% % plot(0.125:0.125:51,beta_a(3,:))
% plot(position,beta(3,:),'-r')
% ylim([0 400])
% xlim([1150 1650])
% 
% testoutput = [position, beta(3,:)', ci(3,:)'];


% plot(position,log10(beta(5,:)),'-r')
%  [wetcoeffs(:,i),wetresnorm(i),wetresidual(:,i)] = lsqcurvefit(@t2monofit_simple_simple,guesses,echoVector,real(wet7decays(:,i)),lb,ub,opts);

%%
%% Do ILT
alpha = 1e8;
k=1;
lowLim = min(final_data(k).echoVector);
hiLim = max(final_data(k).echoVector);
nrILTSteps = length(final_data(k).echoVector);
clear('spectrum','tau','chisq')
spectrum = zeros(nrILTSteps,length(final_data(k).position));
tau = spectrum;
chisq = zeros(length(final_data(k).position),1);
% [spectrum,tau,chisq,~     ] = upnnlsmooth1D(final_data(k).echoVector',final_data(k).data(:,18),  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
% figure
% plot(log10(tau),(real(spectrum)))
%
tic
parfor i = 1:length(final_data(k).position)
%     tic
    [spectrum(:,i),tau(:,i),chisq(i),~     ] = upnnlsmooth1D(final_data(k).echoVector',final_data(k).data(:,i),  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
    disp(strcat('No. ',num2str(i)));
%     toc
end
totalTime = toc;
disp(strcat('Time per loop: ',num2str(totalTime/length(final_data(k).position)),' sec.'));

ILTArea = sum(spectrum,1);
beep;

% rawArea = sum(real(phasedData),1);
% save('processed_data.mat')

% figure
% semilogx(tau,spectrum)

%% Plot ILT
figure(4)
subplot(4,4,[1 2 3 5 6 7 9 10 11])
surf(final_data(1).position,log10(tau),real(spectrum));
view([0 90])
colormap(hot)
% caxis([0 400])
 shading flat
 set(gca,'YDir','Reverse')
 xlabel('position (mm)')
 ylabel('log(T_2/ms)')
 zlabel('intensity')
 ylim([log10(lowLim) log10(hiLim)])
%  xlim([1800 2800])
subplot(4,4,[13 14 15])
plot(final_data(1).position,sum(final_data(1).data));
 set(gca,'YDir','Reverse')
    xlabel('distance (mm) (toward MOUSE <-- --> away from MOUSE)')
%      xlim([1800 2800])
subplot(4,4,[4 8 12])
plot(sum(spectrum'),log10(tau)) %#ok<UDIM>
%  ylim([log10(lowLim) log10(hiLim)])
 line([0 5000],[log10(final_data(1).echoVector(1)/1000) log10(final_data(1).echoVector(1)/1000)])
 line([0 5000],[log10(final_data(1).echoVector(length(final_data(1).echoVector))/1000) log10(final_data(1).echoVector(length(final_data(1).echoVector))/1000)],'LineStyle','-.')
  set(gca,'YDir','Reverse')
%    ylabel('log(T_2/ms)')
