clear
clc
close all

%%

info.parfilename = 'acqu';
info.datafilename = 'Noconsolidant_papyrus';
info.dirstem = '/Volumes/ISC1026/Data/TKM/NGA_12May2016/Hamada/Noconsolidant_papyrus/';
omitpoints = 0;

fitopts = statset('MaxIter',5000,'TolX',1e-14,'UseParallel',true,'Display','off');

dir = strcat(info.dirstem,num2str(3),'/');
cd(dir);

params.acqTime = readpar_Kea(strcat(info.parfilename,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(info.parfilename,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(info.parfilename,'.par'),'nrScans');
params.rxPhase = readpar_Kea(strcat(info.parfilename,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(info.parfilename,'.par'),'rxGain');
params.nrPts = readpar_Kea(strcat(info.parfilename,'.par'),'nrPnts');
params.repTime = readpar_Kea(strcat(info.parfilename,'.par'),'repTime');
params.b1Freq = readpar_Kea(strcat(info.parfilename,'.par'),'b1Freq');
params.nrEchoes = readpar_Kea(strcat(info.parfilename,'.par'),'nrEchoes');
params.echoTime = readpar_Kea(strcat(info.parfilename,'.par'),'echoTime');
params.nrExp = readpar_Kea(strcat(info.parfilename,'.par'),'nrExp');
params.initDepth = readpar_Kea(strcat(info.parfilename,'.par'),'initDepth');
params.finalDepth = readpar_Kea(strcat(info.parfilename,'.par'),'finalDepth');
params.stepSize = readpar_Kea(strcat(info.parfilename,'.par'),'stepSize');

real_data = load(strcat(info.datafilename,'-decaysRe.dat'));
imag_data = load(strcat(info.datafilename,'-decaysIm.dat'));
temp_data = complex(real_data(omitpoints+1:end,2:end),imag_data(omitpoints+1:end,2:end));
data = real(autophase(temp_data,1));
contrast_data = load(strcat(info.datafilename,'.dat'));

%%
final_data.position = contrast_data(:,1);
%final_data.position = 0;
% for i = 1:\length(contrast_data(:,1))
%     final_data.position = vertcat(final_data.position,.4551*(i-1));
% end
final_data.echoVector = params.echoTime*(omitpoints+1):params.echoTime:params.echoTime*params.nrEchoes;
final_data.data = data;

xfit = 0:1:max(final_data(1).echoVector);

guess = [0.07; 340];
%guessesbi = [.5;.5;.4;22];


% Monoexponential fit
for i=1:length(final_data.position)
    
   try
       [final_data.beta(:,i),final_data.residual(:,i),final_data.J(:,:,i)] = nlinfit(final_data.echoVector',final_data.data(:,i),@t2monofit_simple,guess,fitopts);
%                final_data(k).beta(3,i) = 0;
%                final_data(k).beta(4,i) = 0;
   catch pm
       final_data.beta(:,i) = [0;0];%;0;0];
       final_data.J(:,:,i) = zeros(size(final_data.data,1),2);
       final_data.residual(:,i) = zeros(size(final_data.data,1),1);
   end

   %guess(:,i) = final_data.beta(:,i);
   
   [final_data.ypred(:,i),final_data.delta(:,i)] = nlpredci(@t2monofit_simple,xfit,final_data.beta(:,i),final_data.residual(:,i),'Jacobian',final_data.J(:,:,i));
   tt = nlparci(final_data.beta(:,i),final_data.residual(:,i),'jacobian',final_data.J(:,:,i), 'alpha', 0.1);
   final_data.ci(:,i) = tt(:,2)-final_data.beta(:,i);
   [final_data.yvals(:,i),final_data.yvalsdelta(:,i)] = nlpredci(@t2monofit_simple,final_data.echoVector,final_data.beta(:,i),final_data.residual(:,i),'Jacobian',final_data.J(:,:,i));
end

% % Biexponential Fit
% for i=1:length(final_data.position)
%     for ii = 1:10
%        try
%            [final_data.beta(:,i),final_data.residual(:,i),final_data.J(:,:,i)] = nlinfit(final_data.echoVector',final_data.data(:,i),@t2bifit_simple,guessesbi,fitopts);
%            
%        catch
%            final_data.beta(:,i) = [0;0;0;0];
%            final_data.J(:,:,i) = zeros(size(final_data.data,1),4);
%            final_data.residual(:,i) = zeros(size(final_data.data,1),1);
%        end
%        
%        if final_data.beta(:,i) ~= [0;0;0;0]
%            guessesbi = final_data.beta(:,i);
%        else
%            continue
%        end
%        
%     end
%     [final_data.ypred(:,i),final_data.delta(:,i)] = nlpredci(@t2bifit_simple,xfit,final_data.beta(:,i),final_data.residual(:,i),'Jacobian',final_data.J(:,:,i));
%     tt = nlparci(final_data.beta(:,i),final_data.residual(:,i),'jacobian',final_data.J(:,:,i));
%     final_data.ci(:,i) = tt(:,2)-final_data.beta(:,i);
%     [final_data.yvals(:,i),final_data.yvalsdelta(:,i)] = nlpredci(@t2bifit_simple,final_data.echoVector,final_data.beta(:,i),final_data.residual(:,i),'Jacobian',final_data.J(:,:,i));
%     
% end

%%

% for i=1:length(contrast_data(:,1))
%     final_data.position
% 

positionlims = [1 91];
shiftAmount = 80;
figure(1)
subplot(1,2,1)
hold on
plot(final_data.position(positionlims(1):positionlims(2))+shiftAmount,final_data.beta(2,positionlims(1):positionlims(2)),'-o')

% plot(final_data.position(positionlims(1):positionlims(2))+shiftAmount,final_data.beta(2,positionlims(1):positionlims(2))+final_data.ci(2,positionlims(1):positionlims(2)),'--b')
% plot(final_data.position(positionlims(1):positionlims(2))+shiftAmount,final_data.beta(2,positionlims(1):positionlims(2))-final_data.ci(2,positionlims(1):positionlims(2)),'--b')
ylim([-0.05 2000])
xlim([1100 2900])
title('T_2 vs. position')
ylabel('T_2 (us) from fit to subtracted data')
subplot(1,2,2)
hold on
plot(final_data.position(positionlims(1):positionlims(2))+shiftAmount,final_data.beta(1,positionlims(1):positionlims(2)),'-o')
% plot(final_data.position(positionlims(1):positionlims(2))+shiftAmount,final_data.beta(1,positionlims(1):positionlims(2))+final_data.ci(1,positionlims(1):positionlims(2)),'--b')
% plot(final_data.position(positionlims(1):positionlims(2))+shiftAmount,final_data.beta(1,positionlims(1):positionlims(2))-final_data.ci(1,positionlims(1):positionlims(2)),'--b')
 ylim([0 0.3])
xlim([1100 2900])
title('Amplitude vs. position, Dry')
xlabel('canvas <--- position (um) ---> paint')
ylabel('Amplitude (arb) from fit to subtracted data')

% subplot(1,3,3)
% hold on
% plot(final_data.position(positionlims(1):positionlims(2)),RH,'r')
% ylim([20 80])



