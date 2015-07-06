clear
clc
close all

%%


% Load Data and Parameters

%FTdata = load('/Users/jaredking/Documents/Research Files and Data/Spatial Resolution/CPMG/FreshYellow_11Jun/1040_780Series3/980um/1/data2.csv');

% dataRe = load('/Users/jaredking/Documents/Research Files and Data/Spatial Resolution/CPMG/FreshYellow_11Jun/1120_820Series/1000um/1/echoData2D-Re.dat');
% dataIm = load('/Users/jaredking/Documents/Research Files and Data/Spatial Resolution/CPMG/FreshYellow_11Jun/1120_820Series/1000um/1/echoData2D-Im.dat');


parfilestem = sprintf('%s', 'C:\Users\vjlee\Desktop\2\acqu');

params.acqTime = readpar_Kea(strcat(parfilestem,'.par'),'acqTime');
params.bandwidth = readpar_Kea(strcat(parfilestem,'.par'),'bandwidth');
params.nrScans = readpar_Kea(strcat(parfilestem,'.par'),'nrScans');
params.rxPhase = readpar_Kea(strcat(parfilestem,'.par'),'rxPhase');
params.rxGain = readpar_Kea(strcat(parfilestem,'.par'),'rxGain');
params.nrPts = readpar_Kea(strcat(parfilestem,'.par'),'nrPnts');
params.repTime = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.repTime = readpar_Kea(strcat(parfilestem,'.par'),'repTime');
params.b1Freq = readpar_Kea(strcat(parfilestem,'.par'),'b1Freq');
params.nrEchoes = readpar_Kea(strcat(parfilestem,'.par'),'nrEchoes');
params.echoTime = readpar_Kea(strcat(parfilestem,'.par'),'echoTime');

% New Data files
file = sprintf('%s', 'C:\Users\vjlee\Desktop\2\data2.csv');
FTdata = load(file);

%%
FTdata = reshape(FTdata,(params.nrEchoes*2),params.nrPts);
dataRe = FTdata(1:params.nrEchoes,:);
dataIm = FTdata((params.nrEchoes+1):(params.nrEchoes*2),:);

% FTdata = reshape(FTdata,(params.nrEchoes*2),params.nrPts);
% dataRe = FTdata(1:params.nrEchoes,:);
% dataIm = FTdata((params.nrEchoes+1):(params.nrEchoes*2),:);

% %% Fourier Transforms (Tyler's Code)
% 
% % Tweaks final results
% %i = 2; %echo number for plotting
% reduction = 1; %reduce sampling rate by this factor to decrease bandwidth 
% lb = 0; %line broadening
% zerofilling = 2; %factor for ZF (0 = no ZF, 1 = double the length, 2 = 4x the length, etc...)
% sample_window = 50; %window for saving samples in um
% 
% % Goes 'ding' when there's stuff
% T = params.acqTime*(2^zerofilling);     %total sample time (ms)
% L = params.nrPts*(2^zerofilling)/reduction;       %number of complex points
% dt = T/L;               %dwell time (ms)
% Fs = 1/dt;              %sampling freq (kHz)
% t = dt:dt:T;            %time vector
% 
% % Makes space scale
% NFFT = 2^nextpow2(L);
% freq_space = linspace(-Fs/2,Fs/2,NFFT); %kHz
% real_space = freq_space/1016; %mm %Gradient of magnet (852.51)
% 
% g_x = exp(-(t-(max(t)/2)).^2/(2*(1/lb)^2)); % Supresses noise
% Y = zeros(params.nrEchoes,L);
% dY = zeros(params.nrEchoes,L-1);
% 
% % Set complex data
% if zerofilling == 0
%     dataCmpx = complex(dataRe,dataIm);
% else
%     dataCmpx = padarray(complex(dataRe,dataIm), [0 (NFFT-params.nrPts)/2]);
% end
% 
% % % Does the transform for each echo one at a time
% % for j = 1:params.nrEchoes;
% %     data =  sum(dataCmpx(j,:),1); % Changes from complex to abs signal
% %     data_dec = decimate(data,reduction).*g_x; % Decimates data (see above)
% %     Y(j,:) = fftshift(fft(data_dec,NFFT)/L); % Performs the Transform
% %     dY(j,:) = diff(abs(Y(j,:)))./diff(real_space); % Takes derivative of fourier transform
% % end
% 


%% Nlinfit for non-transformed Data

% Sum Normalized Data
for i = 1:params.nrEchoes
    sumData(i,1) = sum(dataRe(i,:));
end
maxVal = max(sumData);
maxVal = max(maxVal);
sumData = sumData/maxVal;

%% Summed Data Monofit

guess = [.9;.65];
for i=1:10
    %                                                       Be sure to change
    %                                                       the echotime
    [beta,Resids,J,covB] = nlinfit((1:params.nrEchoes)'*45e-3,sumData(:,1),@t2monofit_simple,guess);

    % Confidence Interval and Margin of Error (+/-)
    CI=nlparci(beta,Resids,'jacobian',J);
    MOE_T2 = (CI(2,2) - CI(2,1))/2;
    guess = beta;
    
end

%% Summed Data BiFit
echoVector = (1:params.nrEchoes)*45e-3;

guess = [.5;.04;.4;3.07];
for i=1:10
    %                                                       Be sure to change
    %                                                       the echotime
    [beta,Resids,J,covB] = nlinfit(echoVector',sumData,@t2bifit_simple,guess);

    % Confidence Interval and Margin of Error (+/-)
    CI=nlparci(beta,Resids,'jacobian',J);
    MOE_T2 = (CI(2,2) - CI(2,1))/2;
    guess = beta;
    
end

% ypred = t2bifit_simple(beta,echoVector');

% figure(1)
% hold on
% plot(echoVector,sumData);
% plot(echoVector,ypred)

%% Plot Summed Data

figure(7)
plot(sumData(:,:))

%%
% %% T2 Curve Fitting
% 
% % T2 Guess and other Variables
% T2guess = [1;50];
% timeAxis = 1:params.nrEchoes;
% timeAxis = timeAxis' * 160e-3; % change to multiply by pulse length; ms?
% 
% % Calculate mono-exponential fit
% for i=1:NFFT
%     fitdata = abs(Y(:,i));
%     try
%         [beta(:,i), betaResids] = nlinfit(timeAxis,fitdata,@t2monofit_simple,T2guess);
%     catch
%         beta(:,i) = 0:0;
%     end
%         
% end


%%% Fourier Transforms (My Code)

% % Set Complex Data; Plot real and imaginary lines
% if zerofilling == 0
%     dataCmpx = complex(dataRe,dataIm);
% else
%     dataCmpx = padarray(complex(dataRe,dataIm), [0 (NFFT-params.nrPts)/2]);
% end
% 
% %FTfile = fopen('/Users/jaredking/Documents/Research Files and Data/Spatial Resolution/CPMG/HandFTs','a');
% 
% % Perform and plot Fourier Transform
% for i = 1:params.nrEchoes
%     ftpnt = dataCmpx(i,:);
%     Y(i,:) = fftshift(fft(ftpnt));
%     %fprintf(FTfile,'%f\t',abs(Y(i,30:35)));
%     %fprintf(FTfile,'%s\n','');
% end

% for e=1:12
%     for i=1:700
%         fprintf(FTfile,'%s\n','');
%         fprintf(FTfile,'%f\t',abs(Y(i,30:35)));
%     end
% end

%fclose(FTfile);

% %% Sum FT Data and write to file
% FTfile = fopen('/Users/jaredking/Documents/Research Files and Data/Spatial Resolution/CPMG/sumFT_data','a');
% %fread(FTfile);
% sumFT = abs(sum(Y));
% fprintf(FTfile,'%s\n','');
% fprintf(FTfile,'%f\t',sumFT);
% fclose(FTfile);





% %% Plot T2 with Fourier
% 
% figure(6)
% hold on
% %plot(real_space*1000,beta(2,:)/10,'g')
% plot(real_space*1000,abs(Y(1,:)),'r')
% %xlim([-160 80]);
% %ylim([0 10]);
% xlabel('distance (microns)')
% ylabel('time (ms/10)')
% hold off
% 
% %% Plot Single Echo -  Real/Imaginary Data
% figure(1)
% hold on;
% for i = 3
%     plot(real(dataCmpx(i,:)),'b');
% end
% hold on;
% for i = 3
%     plot(imag(dataCmpx(i,:)),'r');
% end
% 
% %% Subplot of each echo
% 
% figure(4)
% hold on;
% for i = 1:params.nrEchoes
%     subplot(4,4,i)
%     hold on;
%     plot(real(dataCmpx(i,:)),'b');
%     plot(imag(dataCmpx(i,:)),'r');
% end
% 
% %% Subplot of Fourier Transforms
% 
% figure(5)
% hold on;
% for i = 1:params.nrEchoes
%     subplot(4,4,i)
%     dataCmpx1 = dataCmpx(i,:);
%     fdCp1 = fftshift(fft(dataCmpx1));  
%     plot(abs(fdCp1),'g');
% end

%%
echoVector = (1:params.nrEchoes)*45e-3;
figure(1)
plot(echoVector',sumData)

ySolv = .5522*exp(-echoVector/.2799)+.5089*exp(-echoVector/4.4285);

figure(2)
hold on
plot(echoVector,ySolv)
plot(echoVector',sumData)
legend('data','solvent')