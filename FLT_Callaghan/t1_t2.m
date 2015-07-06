% for t1- t2 with saturation recovery, t1 changes in log10 scale
% M=Mo*(1-exp{-tau1/T1})*exp{-kte/T2} + noise
% all values should be in secs

clear;

R='L:\Aachen_RWTH\FLT\test\text79';

tau1_st=20;                  % recovery start [ms]
tau1_st=tau1_st/1000;         % -//- [s]

tau1_end=5000.0;                 % recovery end [ms]
tau1_end=tau1_end/1000;       % -//- [s]

te=0.15;                    % echo time [ms]
te=te/1000;


n_T1=164;                     % number of T1
n_echo=128;                  % number of T2
Mo=90;                       % amplitude

T1=1.2;                       % T1 [s]


T2=0.1;                       % T2 [s]

power_incr=log10(tau1_end)-log10(tau1_st);
power_incr=power_incr/n_T1;
incr=10^power_incr;
level_noise=0.0;           % noise level
%noise=rand(n_echo,n_T1); % generate random numbers from 0 till 1
noise=randn(n_echo,n_T1); % generate Gaussian, or normally distributed, random numbers
max_noise=max(abs(noise));
for i=1:n_echo;
    for j=1:n_T1;
        noise(i,j)=noise(i,j)/max_noise(1,j);
    end;
end;

noise=noise*level_noise;

for i=1:n_T1;
    tau1(i,1)=tau1_st*incr^(i-1);
    for j=1:n_echo;
        kte(j,1)=te*j;
        M(j,i)=Mo*1.0*(1-1.0*exp(-tau1(i)/T1))*exp(-te*j/T2);
        M_noise(j,i)=M(j,i);
        M_tr(i,j)=M(j,i);
        M_tr_noise(i,j)=M_noise(j,i);
    end;
end;

Data=M;
Tau_1=tau1;
Tau_2=kte;
NoiseStd=100;

% M=Mmax-M for the callaghan prog
for i=1:n_T1;     
  for j=1:n_echo;
 data_cllghn_I(j,i)=M(j,n_T1)-M(j,i);
  data_cllghn_noise(j,i)=M_noise(j,n_T1)-M_noise(j,i);
   end;
end;

for i=1:n_T1-1;           %  minus posledni nulevaja kolonka
 for j=1:n_echo;     
   Data_minus(j,i)=data_cllghn_I(j,i);
   Data_minus_noise(j,i)=data_cllghn_noise(j,i);
   Data_minus_tr(i,j)=data_cllghn_I(j,i);
   minus_Tau1(i,1)=tau1(i,1);
 end;
end;
 
plot(tau1,M);figure(gcf)
%figure1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
savefile=[R,'_sat.jpg'];
saveas(gcf,savefile);
clear savefile;

plot(tau1,M_noise);figure(gcf)
%figure1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
savefile=[R,'_sat_noise.jpg'];
saveas(gcf,savefile);
clear savefile;


plot(kte,M_tr);figure(gcf)
%figure1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
savefile=[R,'_cpmg.jpg'];
saveas(gcf,savefile);
clear savefile;

plot(kte,M_tr_noise);figure(gcf)
savefile=[R,'_cpmg_noise.jpg'];
saveas(gcf,savefile);
clear savefile;

clear i;
clear j;
clear power_incr;
clear figure1;

%SaveFile=strcat(R,'_simulation_data.mat');
%save (SaveFile);
%clear SaveFile;

% for callaghan
%save F:\simulations\ac6\t1t2\te60us\ac6_minus_data Data_minus -ASCII;

SaveFile=strcat(R,'_minus_data_s');
%dlmwrite (SaveFile, Data_minus_noise, ' ');
save ( SaveFile, 'Data_minus_noise','-ASCII')
clear SaveFile;
%fprintf('now amplitudes........................now amplitudes ');

SaveFile=strcat(R,'_tau2');
csvwrite (SaveFile , kte);
clear SaveFile;
%fprintf('Echoes are saved :D:D:D:D');

SaveFile=strcat(R,'_tau1_s');
csvwrite (SaveFile, minus_Tau1);
clear SaveFile;


% for huerlimann 
SaveFile=strcat(R,'_Huerlm.mat');
save (SaveFile, 'Tau_1',  'Tau_2',  'Data', 'NoiseStd');
clear SaveFile;

fprintf('FINISHED!!!!......FINISHED!!!........FINISHED!!!.......FINISHED!!! ');




