clear
clc
close all
% try
%     matlabpool
% catch
% end
% Load data
filedir = '/Users/tyler/Dropbox/Data/NGA/NGA_5June2015/Atlas9pARLO_1024/1/';

cd(filedir)
filename = 'data.csv';
parname = 'acqu';
alpha = 1e10;
omitpoints = 0;

data = load(filename);
echoVector = data(omitpoints+1:end,1)/1e6;
% real_data =  data(omitpoints+1:end,2);
% imag_data =  data(omitpoints+1:end,3);
realData = data(omitpoints+1:end,2)./max(data(omitpoints+1:end,2));

% real_data = load(strcat(filename,'-decaysRe.dat'));
% imag_data = load(strcat(filename,'-decaysIm.dat'));
% 
% contrast_data = load(strcat(filename,'.dat'));
% position = contrast_data(:,1);

% cplx_data = complex(real_data,imag_data);
% phased_data = autophase(cplx_data,0.1);

% params.acqTime = readpar_Kea(strcat(parname,'.par'),'acqTime');
% params.bandwidth = readpar_Kea(strcat(parname,'.par'),'bandwidth');
% params.nrScans = readpar_Kea(strcat(parname,'.par'),'nrScansT2');
% params.rxPhase = readpar_Kea(strcat(parname,'.par'),'rxPhase');
% params.rxGain = readpar_Kea(strcat(parname,'.par'),'rxGain');
% params.nrPts = readpar_Kea(strcat(parname,'.par'),'nrPnts');
% params.repTime = readpar_Kea(strcat(parname,'.par'),'repTimeT2');
% params.b1Freq = readpar_Kea(strcat(parname,'.par'),'b1Freq');
% params.nrEchoes = readpar_Kea(strcat(parname,'.par'),'nrEchoesT2');
% params.echoTime = readpar_Kea(strcat(parname,'.par'),'echoTime');

% echoVector = params.echoTime:params.echoTime:params.echoTime*params.nrEchoes;
lowLim = 10^-5; %min(echoVector)/10000; %
hiLim = 10^-2; %max(echoVector)/10;
nrILTSteps = length(echoVector)/1;

scatter(echoVector,realData)

[spectrum,tau,chisq,~     ] = upnnlsmooth1D(realData,echoVector,  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
semilogx(tau,spectrum)
xlabel('T_2 time [s]')
ylabel('intensity [arb]')
ILTout = [tau',spectrum'];

save('ILTout.dat','ILTout','-ascii');

%% DO a 1D plot and ILT on the strongest signal to test for multiexponentiality
% [C,I] = max(real(cplx_data(1,:)));
% 
% guesses = [ 0;
%             20;
%             1;
%             150;
%             0.01];
% CI = 90;
% 
% [spectrum1D,tau1D,chisq1D,~     ] = upnnlsmooth1D(real(cplx_data(:,I)),echoVector'/1000,  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
% 
% [xfit1D,ypred1D,beta1D,betaErr1D,resid1D] = bidecay_t2fit(echoVector'/1000,real(cplx_data(:,I)),guesses,CI);
% 
% figure(1)
% subplot(2,1,1)
% hold on
% plot(echoVector/1000,real(cplx_data(:,I)),'-k')
% plot(echoVector/1000,imag(cplx_data(:,I)),'-r')
% plot(xfit1D,ypred1D,'-b')
% xlabel('time (ms)')
% xlim([min(echoVector/1000) max(echoVector/1000)])
% subplot(2,1,2)
% hold on
% plot(log10(tau1D),spectrum1D)
% line([log10(beta1D(3,1)) log10(beta1D(3,1))],[0 30],'LineStyle','-')
% line([log10(beta1D(3,1)-betaErr1D(3)) log10(beta1D(3,1)-betaErr1D(3))],[0 30],'LineStyle','--')
% line([log10(beta1D(3,1)+betaErr1D(3)) log10(beta1D(3,1)+betaErr1D(3))],[0 30],'LineStyle','--')
% line([log10(beta1D(5,1)) log10(beta1D(5,1))],[0 30],'LineStyle','-','Color','r')
% line([log10(beta1D(5,1)-betaErr1D(5)) log10(beta1D(5,1)-betaErr1D(5))],[0 30],'LineStyle','--','Color','r')
% line([log10(beta1D(5,1)+betaErr1D(5)) log10(beta1D(5,1)+betaErr1D(5))],[0 30],'LineStyle','--','Color','r')
% xlabel('log time/ms')
% xlim([log10(min(echoVector/1000)) log10(hiLim)])
% 
% %%
% parfor i = 1:length(position)
%     try
%     [xfit(:,i),ypred(:,i),beta(:,i),beta_err(:,i),resid(:,i)] = bidecay_t2fit(echoVector'/1000,real(cplx_data(:,i)),guesses,CI);
%     catch 
% %         pm
% %         beta(:,i) = 0;
% %         beta_err(:,i) = 0;
% %         resid(:,i) = 0;
%     end
% end

%%
% clear('spectrum','tau','chisq')
% spectrum = zeros(nrILTSteps,1);
% tau = spectrum;
% chisq = zeros(length(position));

% tic

%     parfor i = 1:length(position)
        %     tic
        [spectrum,tau,chisq,~     ] = upnnlsmooth1D(real(cplx_data),echoVector,  lowLim, hiLim, alpha ,  -1,  nrILTSteps);
%         disp(strcat('No. ',num2str(i)));
        %     toc
%     end

semilogx(tau,spectrum)
% totalTime = toc;
% disp(strcat('Time per loop: ',num2str(totalTime/length(position)),' sec.'));

% %% Plot
% 
% beta_relErr1 = beta_err(3,:)./beta(3,:);
% beta_select1 = (abs(beta_relErr1)<1);
% beta_relErr2 = beta_err(5,:)./beta(5,:);
% beta_select2 = (abs(beta_relErr2)<1);
% beta(3,:) = beta_select1.*beta(3,:);
% beta_err(3,:) = beta_err(3,:).*beta_select1;
% beta(5,:) = beta_select2.*beta(5,:);
% beta_err(5,:) = beta_err(5,:).*beta_select2;
% 
% p = figure(4);
% subplot(4,4,[1 2 3 5 6 7 9 10 11])
% title('Untreated Insert 2')
% hold on
% surf(position,log10(tau(:,1,1)),spectrum)
% % surf(position,tau(:,1,1),spectrum)
% shading interp
% % contourf(position,log10(tau(:,1,1)),spectrum,[5 10 15 20 25]);
% % h = findobj('Type','patch');
% % set(h,'LineWidth',1)
% % plot(position,log10(beta(3,:)),'-b');
% % plot(position,log10(beta(3,:)-beta_err(3,:)),':b');
% % plot(position,log10(beta(3,:)+beta_err(3,:)),':b');
% % plot(position,log10(beta(5,:)),'-r')
% % plot(position,log10(beta(5,:)-beta_err(5,:)),':r');
% % plot(position,log10(beta(5,:)+beta_err(5,:)),':r');
% view([0 90])
% revgray = [linspace(1,0)', linspace(1,0)', linspace(1,0)'];
% colormap(revgray)
% caxis([0 20])
% %  shading flat
%  set(gca,'YDir','Reverse')
%  xlabel('position (mm)')
%  ylabel('log(T_2/ms)')
%  zlabel('intensity')
% %  ylim([log10(echoVector(1)/1000) log10(hiLim)])
% ylim([-1.5 0.5])
%  zlim([0 140])
%  xlim([1800 4000])
% subplot(4,4,[13 14 15])
% plot(position,sum(real(cplx_data)));
%  set(gca,'YDir','Reverse')
%     xlabel('distance (mm) (toward MOUSE <-- --> away from MOUSE)')
%      xlim([1800 4000])
% subplot(4,4,[4 8 12])
% plot(sum(spectrum'),log10(tau(:,:,1)))
% % plot(sum(spectrum'),tau(:,:,1))
% %  ylim([log10(echoVector(1)/1000) log10(hiLim)])
% %  line([0 5000],[log10(echoVector(1)/1000) log10(echoVector(1)/1000)])
%  line([0 5000],[log10(echoVector(length(echoVector))/1000) log10(echoVector(length(echoVector))/1000)],'LineStyle','-.')
%   set(gca,'YDir','Reverse')
%   ylim([-1.5 0.5])
%   xlim([0 2000])
% %    ylabel('log(T_2/ms)')
%%
% orient landscape
% print(p,'-dpdf','UntreatedInsert','-r300')
% print(p,'-depsc2','UntreatedInsert','-r300')
% print(p,'-dtiff','UntreatedInsert','-r300')
